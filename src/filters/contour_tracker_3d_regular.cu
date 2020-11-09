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

using namespace ftk;

template <int scope>
__device__ __host__
bool check_simplex_contour_3dt(
    int current_timestep,
    const lattice4_t& domain, 
    const lattice4_t& core, 
    const lattice3_t& ext, // array dimension
    const element41_t& e,
    double threshold,
    const double *F[2],
    cp_t &p)
{
  if (e.corner[3] != current_timestep)
    return false;
  
  int vertices[2][4], indices[2];
  size_t local_indices[2];
  for (int i = 0; i < 2; i ++) {
    for (int j = 0; j < 4; j ++) {
      vertices[i][j] = e.corner[j]
        + unit_simplex_offset_4_1<scope>(e.type, i, j);
      if (vertices[i][j] < domain.st[j] || 
          vertices[i][j] > domain.st[j] + domain.sz[j] - 1)
        return false;
    }
    indices[i] = domain.to_index(vertices[i]);
    local_indices[i] = ext.to_index(vertices[i]);
  }
  
  const long long factor = 2 << 20; // using fixed point rep
  double X[2][4];
  double f[2];
  long long fi[2];
  for (int i = 0; i < 2; i ++) {
    const size_t k = local_indices[i]; // k = ext.to_index(vertices[i]);
    const size_t t = unit_simplex_offset_4_1<scope>(e.type, i, 3);
     
    f[i] = F[t][k] - threshold;
    fi[i] = f[i] * factor;

    for (int j = 0; j < 4; j ++)
      X[i][j] = vertices[i][j];
  }
 
  bool succ = robust_critical_point_in_simplex1(fi, indices);
  if (!succ) return false;

  double mu[2];
  bool succ2 = inverse_lerp_s1v1(f, mu);
  
  double x[4];
  lerp_s1v4(X, mu, x);

  p.x[0] = x[0];
  p.x[1] = x[1];
  p.x[2] = x[2];
  p.t = x[3];

  return true;
}

template <int scope>
__global__
void sweep_simplices(
    int current_timestep,
    const lattice4_t domain,
    const lattice4_t core,
    const lattice3_t ext, // array dimension
    double threshold,
    const double *F_c, // current timestep
    const double *F_n, // next timestep
    unsigned long long &ncps, cp_t *cps)
{
  const double *F[2] = {F_c, F_n};
  
  int tid = getGlobalIdx_3D_1D();
  const element41_t e = element41_from_index<scope>(core, tid);
  
  cp_t cp;
  bool succ = check_simplex_contour_3dt<scope>(
      current_timestep,
      domain, core, ext, e, threshold, F, cp);

  if (succ) {
    unsigned long long i = atomicAdd(&ncps, 1ul);
    cp.tag = tid;
    cps[i] = cp;
  }
}

template <int scope>
static std::vector<cp_t> extract_contour_3dt(
    int current_timestep,
    const lattice4_t& domain,
    const lattice4_t& core, 
    const lattice3_t& ext,
    double threshold,
    const double *F_c,
    const double *F_n)
{
  auto t0 = std::chrono::high_resolution_clock::now();

  const size_t ntasks = core.n() * ntypes_4_1<scope>();
  // fprintf(stderr, "ntasks=%zu\n", ntasks);
  const int maxGridDim = 1024;
  const int blockSize = 256;
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;
  
  if (nBlocks >= maxGridDim) 
    gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else 
    gridSize = dim3(nBlocks);

  double *dF_c = NULL, *dF_n = NULL;
  if (F_c) {
    cudaMalloc((void**)&dF_c, sizeof(double) * ext.n());
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: allocating dF_c");
    cudaMemcpy(dF_c, F_c, sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: copying drho_c");
  }
  if (F_n) {
    cudaMalloc((void**)&dF_n, sizeof(double) * ext.n());
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: allocating drho_l");
    cudaMemcpy(dF_n, F_n, sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: copying drho_l");
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
      domain, core, ext, threshold, dF_c, dF_n, 
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
 
  if (dF_c) cudaFree(dF_c);
  if (dF_n) cudaFree(dF_n);
  cudaFree(dncps);
  cudaFree(dcps);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaFree");
 
  cudaDeviceSynchronize();
  auto t1 = std::chrono::high_resolution_clock::now();
  double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() * 1e-9;

  fprintf(stderr, "exitting gpu kernel, ncps=%llu, time=%f\n", ncps, duration);
  
  return cps;
}

std::vector<cp_t>
extract_contour_3dt_cuda(
    int scope, 
    int current_timestep, 
    const ftk::lattice& domain,
    const ftk::lattice& core, 
    const ftk::lattice& ext, 
    double threshold,
    const double *F_c,
    const double *F_n)
{
  lattice4_t D(domain);
  lattice4_t C(core);
  lattice3_t E(ext);

  if (scope == scope_interval) 
    return extract_contour_3dt<scope_interval>(current_timestep, D, C, E, threshold, F_c, F_n);
  else
    return extract_contour_3dt<scope_ordinal>(current_timestep, D, C, E, threshold, F_c, F_n);
}
