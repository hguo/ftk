#include <nvfunctional>
#include <cstdio>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/symmetric_matrix.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/hypermesh/lattice.hh>
#include <ftk/filters/critical_point.hh>
#include "common.cuh"

//// 
template <int scope>
__device__
bool check_simplex_cp2t(
    int current_timestep,
    const lattice3_t& domain,
    const lattice3_t& core, 
    const lattice2_t& ext, 
    const element32_t& e, 
    const double *V[2], // last and current timesteps
    const double *gradV[2], // jacobian of last and current timesteps
    const double *scalar[2],
    cp3_t &cp)
{
  typedef ftk::fixed_point<> fp_t;

  const int last_timestep = current_timestep - 1;
  if (scope == scope_interval && e.corner[2] != last_timestep)
    return false;

  int vertices[3][3], indices[3];
  for (int i = 0; i < 3; i ++) {
    for (int j = 0; j < 3; j ++) {
      vertices[i][j] = e.corner[j] 
        + unit_simplex_offset_3_2<scope>(e.type, i, j);
      if (vertices[i][j] < domain.st[j] || 
          vertices[i][j] > domain.st[j] + domain.sz[j] - 1)
        return false;
    }
    indices[i] = domain.to_index(vertices[i]);;
  }

  double v[3][2];
  fp_t vf[3][2];
  for (int i = 0; i < 3; i ++) {
    size_t k = ext.to_index(vertices[i]);
    for (int j = 0; j < 2; j ++) {
      v[i][j] = V[unit_simplex_offset_3_2<scope>(e.type, i, 2/*time dimension id*/)][k*2+j];
      vf[i][j] = v[i][j];
    }
  }

  double mu[3];
  bool succ = ftk::inverse_lerp_s2v2(v, mu, 0.0);
 
  if (succ) {
    // linear jacobian interpolation
    if (gradV[1]) { // have given jacobian
      double Js[3][2][2], J[2][2];
      for (int i = 0; i < 3; i ++) {
        size_t ii = ext.to_index(vertices[i]);
        int t = unit_simplex_offset_3_2<scope>(e.type, i, 2);
        for (int j = 0; j < 2; j ++) 
          for (int k = 0; k < 2; k ++)
            Js[i][j][k] = gradV[t][ii*4 + j*2 + k];
      }
      ftk::lerp_s2m2x2(Js, mu, J);
      ftk::make_symmetric2x2(J);
      cp.type = ftk::critical_point_type_2d(J, true/*symmetric*/);
      // if (cp.type != 0x100) return false;
    }

    // scalar interpolation
    if (scalar[1]) { // have given scalar
      double values[3];
      for (int i = 0; i < 3; i ++) {
        size_t ii = ext.to_index(vertices[i]);
        int t = unit_simplex_offset_3_2<scope>(e.type, i, 2);
        values[i] = scalar[t][ii];
      }
      cp.scalar = ftk::lerp_s2(values, mu);
      // if (abs(cp.scalar) < 0.02) return false; // threshold
    }

    // location interpolation
    double X[3][3];
    for (int i = 0; i < 3; i ++)
      for (int j = 0; j < 3; j ++)
        X[i][j] = vertices[i][j];
    ftk::lerp_s2v3(X, mu, cp.x);

    return true;
  } else 
    return false;
}

template <int scope>
__global__
void sweep_simplices(
    int current_timestep,
    const lattice3_t domain,
    const lattice3_t core,
    const lattice2_t ext, // array dimensions
    const double *Vc, // current timestep
    const double *Vl, // last timestep
    const double *Jc, 
    const double *Jl,
    const double *Sc, 
    const double *Sl,
    unsigned long long &ncps, cp3_t *cps)
{
  const double *V[2] = {Vl, Vc};
  const double *J[2] = {Jl, Jc};
  const double *S[2] = {Sl, Sc};
  
  int tid = getGlobalIdx_3D_1D();
  const element32_t e = element32_from_index<scope>(core, tid);

  cp3_t cp;
  bool succ = check_simplex_cp2t<scope>(
      current_timestep, 
      domain, core, ext, e, V, J, S, cp);

  if (succ) {
    unsigned long long i = atomicAdd(&ncps, 1ul);
    cp.tag = tid;
    cps[i] = cp;
  }
}

template<int scope>
static std::vector<cp3_t> extract_cp2dt(
    int current_timestep,
    const lattice3_t& domain,
    const lattice3_t& core, 
    const lattice2_t& ext, 
    const double *Vc, // 3D array: 2*W*H
    const double *Vl, 
    const double *Jc,
    const double *Jl,
    const double *Sc,
    const double *Sl)
{
  const size_t ntasks = core.n() * ntypes_3_2<scope>();
  const int maxGridDim = 1024;
  const int blockSize = 256;
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;
  
  if (nBlocks >= maxGridDim) 
    gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else 
    gridSize = dim3(nBlocks);

  double *dVc, *dVl = NULL;
  cudaMalloc((void**)&dVc, 2 * sizeof(double) * ext.n());
  cudaMemcpy(dVc, Vc, 2 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&dVl, 2 * sizeof(double) * ext.n());
  cudaMemcpy(dVl, Vl, 2 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);

  double *dJc = NULL, *dJl = NULL;
  if (Jc) {
    cudaMalloc((void**)&dJc, 4 * sizeof(double) * ext.n());
    cudaMemcpy(dJc, Jc, 4 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }
  if (Jl) {
    cudaMalloc((void**)&dJl, 4 * sizeof(double) * ext.n());
    cudaMemcpy(dJl, Jl, 4 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
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

  cp3_t *dcps;
  cudaMalloc((void**)&dcps, sizeof(cp3_t) * ext.n());
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMalloc/cudaMemcpy");

  fprintf(stderr, "calling kernel func...\n");
  sweep_simplices<scope><<<gridSize, blockSize>>>(
      current_timestep, 
      domain, core, ext, dVc, dVl, dJc, dJl, dSc, dSl,
      *dncps, dcps);
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices");

  unsigned long long ncps;
  cudaMemcpy(&ncps, dncps, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMemcpy, ncps");
  fprintf(stderr, "ncps=%lu\n", ncps);

  std::vector<cp3_t> cps(ncps);
  cudaMemcpy(cps.data(), dcps, sizeof(cp3_t) * ncps, cudaMemcpyDeviceToHost);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMemcpy, dcps");
  
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
  fprintf(stderr, "exitting gpu kernel, ncps=%lu\n", ncps);

  return cps;
}

std::vector<cp3_t>
extract_cp2dt_cuda(
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
  lattice3_t D(domain);
  lattice3_t C(core);
  lattice2_t E(ext);

  // std::cerr << "domain=" << domain 
  //   << ", core=" << core << ", current_timestep=" 
  //   << current_timestep << std::endl;

  if (scope == scope_interval) 
    return extract_cp2dt<scope_interval>(current_timestep, D, C, E, Vc, Vl, Jc, Jl, Sc, Sl);
  if (scope == scope_ordinal) 
    return extract_cp2dt<scope_ordinal>(current_timestep, D, C, E, Vc, Vl, Jc, Jl, Sc, Sl);
  else // scope == 2
    return extract_cp2dt<scope_all>(current_timestep, D, C, E, Vc, Vl, Jc, Jl, Sc, Sl);
}
