#include <nvfunctional>
#include <cstdio>
#include <cassert>
#include <chrono>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/clamp.hh>
#include <ftk/numeric/symmetric_matrix.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/mesh/lattice.hh>
// #include <ftk/filters/critical_point_lite.hh>
#include "common.cuh"

template <int scope>
__device__
bool check_simplex_cp3t(
    int current_timestep,
    const lattice4_t& domain, 
    const lattice4_t& core, 
    const lattice3_t& ext, // array dimension
    const element43_t& e, 
    const double *V[2], // current and next timesteps
    const double *gradV[2], // jacobians
    const double *scalar[2], // scalars
    cp_t &cp)
{
  typedef ftk::fixed_point<> fp_t;
  
  // const int last_timestep = current_timestep - 1;
  // if (scope == scope_interval && e.corner[3] != last_timestep)
  if (e.corner[3] != current_timestep)
    return false;
  
  int vertices[4][4], indices[4];
  size_t local_indices[4];
  for (int i = 0; i < 4; i ++) {
    for (int j = 0; j < 4; j ++) {
      vertices[i][j] = e.corner[j]
        + unit_simplex_offset_4_3<scope>(e.type, i, j);
      if (vertices[i][j] < domain.st[j] || 
          vertices[i][j] > domain.st[j] + domain.sz[j] - 1)
        return false;
    }
    indices[i] = domain.to_index(vertices[i]);
    local_indices[i] = ext.to_index(vertices[i]);
  }

  double v[4][3];
  fp_t vf[4][3];
  for (int i = 0; i < 4; i ++) {
    const size_t k = local_indices[i]; // k = ext.to_index(vertices[i]);
    for (int j = 0; j < 3; j ++) {
      v[i][j] = V[unit_simplex_offset_4_3<scope>(e.type, i, 3)][k*3+j]; // V has three channels
      vf[i][j] = v[i][j];
    }
  }

  bool succ = robust_critical_point_in_simplex3(vf, indices);
  if (!succ) return false;

  double mu[4], cond;
  bool succ2 = ftk::inverse_lerp_s3v3(v, mu, &cond); //, 0.0);
  if (!succ2) ftk::clamp_barycentric<4>(mu);

  if (1) { // (succ2) {
    // if (!succ2) return false;
    // linear jacobian interpolation
    if (gradV[0]) { // have given jacobian
      double Js[4][3][3], J[3][3];
      for (int i = 0; i < 4; i ++) {
        size_t ii = local_indices[i]; // ext.to_index(vertices[i]);
        int t = unit_simplex_offset_4_3<scope>(e.type, i, 3);
        for (int j = 0; j < 3; j ++) 
          for (int k = 0; k < 3; k ++)
            Js[i][j][k] = gradV[t][ii*9 + j*3 + k];
      }
      ftk::lerp_s3m3x3(Js, mu, J);
      ftk::make_symmetric3x3(J);
      cp.type = ftk::critical_point_type_3d(J, true/*symmetric*/);
    }

    // scalar interpolation
    if (scalar[0]) { // have given scalar
      double values[4];
      for (int i = 0; i < 4; i ++) {
        size_t ii = local_indices[i]; // ext.to_index(vertices[i]);
        int t = unit_simplex_offset_4_3<scope>(e.type, i, 3);
        values[i] = scalar[t][ii];
      }
      cp.scalar[0] = ftk::lerp_s3(values, mu);
    }

    double X[4][4], x[4];
    for (int i = 0; i < 4; i ++)
      for (int j = 0; j < 4; j ++)
        X[i][j] = vertices[i][j];
    ftk::lerp_s3v4(X, mu, x);
    cp.x[0] = x[0];
    cp.x[1] = x[1];
    cp.x[2] = x[2];
    cp.t = x[3];
    // cp.cond = cond;
    return true;
  } else 
    return false;
}

template <int scope>
__global__
void sweep_simplices(
    int current_timestep,
    const lattice4_t domain,
    const lattice4_t core,
    const lattice3_t ext, // array dimension
    const double *Vc, // current timestep
    const double *Vn, // next timestep
    const double *Jc, 
    const double *Jn,
    const double *Sc, 
    const double *Sn,
    unsigned long long &ncps, cp_t *cps)
{
  const double *V[2] = {Vc, Vn};
  const double *J[2] = {Jc, Jn};
  const double *S[2] = {Sc, Sn};
  
  int tid = getGlobalIdx_3D_1D();
  const element43_t e = element43_from_index<scope>(core, tid);

  cp_t cp;
  bool succ = check_simplex_cp3t<scope>(
      current_timestep,
      domain, core, ext, e, V, J, S, cp);

  if (succ) {
    unsigned long long i = atomicAdd(&ncps, 1ul);
    cp.tag = tid;
    cps[i] = cp;
  }
}

template <int scope>
static std::vector<cp_t> extract_cp3dt(
    int current_timestep,
    const lattice4_t& domain,
    const lattice4_t& core, 
    const lattice3_t& ext, 
    const double *Vc, 
    const double *Vn, 
    const double *Jc,
    const double *Jn,
    const double *Sc,
    const double *Sn)
{
  auto t0 = std::chrono::high_resolution_clock::now();

  const size_t ntasks = core.n() * ntypes_4_3<scope>();
  // fprintf(stderr, "ntasks=%zu\n", ntasks);
  const int maxGridDim = 1024;
  const int blockSize = 256;
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;
  
  if (nBlocks >= maxGridDim) 
    gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else 
    gridSize = dim3(nBlocks);

  double *dVc = NULL, *dVn = NULL;
  if (Vc) {
    cudaMalloc((void**)&dVc, 3 * sizeof(double) * ext.n());
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: allocating dVc");
    cudaMemcpy(dVc, Vc, 3 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: copying dVc");
  }
  if (Vn) {
    cudaMalloc((void**)&dVn, 3 * sizeof(double) * ext.n());
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: allocating dVl");
    cudaMemcpy(dVn, Vn, 3 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: copying dVl");
  }
  
  double *dJc = NULL, *dJn = NULL;
  if (Jc) {
    cudaMalloc((void**)&dJc, 9 * sizeof(double) * ext.n());
    cudaMemcpy(dJc, Jc, 9 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }
  if (Jn) {
    cudaMalloc((void**)&dJn, 9 * sizeof(double) * ext.n());
    cudaMemcpy(dJn, Jn, 9 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }
  
  double *dSc = NULL, *dSn = NULL;
  if (Sc) {
    cudaMalloc((void**)&dSc, sizeof(double) * ext.n());
    cudaMemcpy(dSc, Sc, sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }
  if (Sn) {
    cudaMalloc((void**)&dSn, sizeof(double) * ext.n());
    cudaMemcpy(dSn, Sn, sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }

  unsigned long long *dncps; // number of cps
  cudaMalloc((void**)&dncps, sizeof(unsigned long long));
  cudaMemset(dncps, 0, sizeof(unsigned long long));
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: allocating dncps");

  cp_t *dcps;
  cudaMalloc((void**)&dcps, sizeof(cp_t) * core.n());
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: allocating dcps");
  cudaDeviceSynchronize();

  fprintf(stderr, "calling kernel func...\n");
  sweep_simplices<scope><<<gridSize, blockSize>>>(
      current_timestep, 
      domain, core, ext, dVc, dVn, dJc, dJn, dSc, dSn,
      *dncps, dcps);
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices, kernel function");

  unsigned long long ncps = 0;
  cudaMemcpy(&ncps, dncps, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMemcpyDeviceToHost, dncps");
  fprintf(stderr, "ncps=%llu\n", ncps);

  std::vector<cp_t> cps(ncps);
  cudaMemcpy(cps.data(), dcps, sizeof(cp_t) * ncps, cudaMemcpyDeviceToHost);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMemcpyDeviceToHost");
  
  if (dVc) cudaFree(dVc);
  if (dVn) cudaFree(dVn);
  if (dJc) cudaFree(dJc);
  if (dJn) cudaFree(dJn);
  if (dSc) cudaFree(dSc);
  if (dSn) cudaFree(dSn);
  cudaFree(dncps);
  cudaFree(dcps);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaFree");
 
  cudaDeviceSynchronize();
  auto t1 = std::chrono::high_resolution_clock::now();
  float duration = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() * 1e-9;

  fprintf(stderr, "exitting gpu kernel, ncps=%llu, time=%f\n", ncps, duration);
  

  return cps;
}

std::vector<cp_t>
extract_cp3dt_cuda(
    int scope, 
    int current_timestep, 
    const ftk::lattice& domain,
    const ftk::lattice& core, 
    const ftk::lattice& ext, 
    const double *Vc, 
    const double *Vl,
    const double *Jc, 
    const double *Jl, 
    const double *Sc,
    const double *Sl)
{
  lattice4_t D(domain);
  lattice4_t C(core);
  lattice3_t E(ext);

  if (scope == scope_interval) 
    return extract_cp3dt<scope_interval>(current_timestep, D, C, E, Vc, Vl, Jc, Jl, Sc, Sl);
  if (scope == scope_ordinal) 
    return extract_cp3dt<scope_ordinal>(current_timestep, D, C, E, Vc, Vl, Jc, Jl, Sc, Sl);
  else // scope == 2
    return extract_cp3dt<scope_all>(current_timestep, D, C, E, Vc, Vl, Jc, Jl, Sc, Sl);
}
