#include <nvfunctional>
#include <cstdio>
#include <cassert>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/symmetric_matrix.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/mesh/lattice.hh>
#include <ftk/filters/critical_point.hh>
#include "common.cuh"

template <int scope>
__device__
bool check_simplex_cp3t(
    int current_timestep,
    const lattice4_t& domain, 
    const lattice4_t& core, 
    const lattice3_t& ext, // array dimension
    const element43_t& e, 
    const double *V[2], // last and current timesteps
    const double *gradV[2], // jacobian of last and current timesteps
    const double *scalar[2],
    cp_t &cp)
{
  const int last_timestep = current_timestep - 1;
  if (scope == scope_interval && e.corner[3] != last_timestep)
    return false;
  
  int vertices[4][4];
  for (int i = 0; i < 4; i ++)
    for (int j = 0; j < 4; j ++) {
      vertices[i][j] = e.corner[j]
        + unit_simplex_offset_4_3<scope>(e.type, i, j);
      if (vertices[i][j] < domain.st[j] || 
          vertices[i][j] > domain.st[j] + domain.sz[j] - 1)
        return false;
    }

  double v[4][3];
  for (int i = 0; i < 4; i ++) {
    size_t k = ext.to_index(vertices[i]);
    for (int j = 0; j < 3; j ++)
      v[i][j] = V[unit_simplex_offset_4_3<scope>(e.type, i, 3)][k*3+j]; // V has three channels
  }

  double mu[4];
  bool succ = ftk::inverse_lerp_s3v3(v, mu);
 
  if (succ) {
    // linear jacobian interpolation
    if (gradV[1]) { // have given jacobian
      double Js[4][3][3], J[3][3];
      for (int i = 0; i < 4; i ++) {
        size_t ii = ext.to_index(vertices[i]);
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
    if (scalar[1]) { // have given scalar
      double values[4];
      for (int i = 0; i < 4; i ++) {
        size_t ii = ext.to_index(vertices[i]);
        int t = unit_simplex_offset_3_2<scope>(e.type, i, 3);
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
    const double *Vl, // last timestep
    const double *Jc, 
    const double *Jl,
    const double *Sc, 
    const double *Sl,
    unsigned long long &ncps, cp_t *cps)
{
  const double *V[2] = {Vl, Vc};
  const double *J[2] = {Jl, Jc};
  const double *S[2] = {Sl, Sc};
  
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
    const double *Vl, 
    const double *Jc,
    const double *Jl,
    const double *Sc,
    const double *Sl)
{
  const size_t ntasks = core.n() * ntypes_4_3<scope>();
  const int maxGridDim = 1024;
  const int blockSize = 256;
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;
  
  if (nBlocks >= maxGridDim) 
    gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else 
    gridSize = dim3(nBlocks);

  double *dVc, *dVl;
  cudaMalloc((void**)&dVc, 3 * sizeof(double) * ext.n());
  cudaMemcpy(dVc, Vc, 3 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&dVl, 3 * sizeof(double) * ext.n());
  cudaMemcpy(dVl, Vl, 3 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  
  double *dJc = NULL, *dJl = NULL;
  if (Jc) {
    cudaMalloc((void**)&dJc, 9 * sizeof(double) * ext.n());
    cudaMemcpy(dJc, Jc, 9 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }
  if (Jl) {
    cudaMalloc((void**)&dJl, 9 * sizeof(double) * ext.n());
    cudaMemcpy(dJl, Jl, 9 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }
  
  double *dSc = NULL, *dSl = NULL;
  if (Sc) {
    cudaMalloc((void**)&dSc, sizeof(double) * ext.n());
    cudaMemcpy(dSc, Sc, sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }
  if (Sl) {
    cudaMalloc((void**)&dSl, sizeof(double) * ext.n());
    cudaMemcpy(dSl, Sl, sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }

  unsigned long long *dncps; // number of cps
  cudaMalloc((void**)&dncps, sizeof(unsigned long long));
  cudaMemset(dncps, 0, sizeof(unsigned long long));

  cp_t *dcps;
  cudaMalloc((void**)&dcps, sizeof(cp_t) * core.n());
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMalloc/cudaMemcpy");

  fprintf(stderr, "calling kernel func...\n");
  sweep_simplices<scope><<<gridSize, blockSize>>>(
      current_timestep, 
      domain, core, ext, dVc, dVl, dJc, dJl, dSc, dSl,
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
  
  cudaFree(dVc);
  cudaFree(dVl);
  if (dJc) cudaFree(dJc);
  if (dJl) cudaFree(dJl);
  if (dSc) cudaFree(dSc);
  if (dSl) cudaFree(dSl);
  cudaFree(dncps);
  cudaFree(dcps);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaFree");
 
  cudaDeviceSynchronize();
  fprintf(stderr, "exitting gpu kernel, ncps=%llu\n", ncps);

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
