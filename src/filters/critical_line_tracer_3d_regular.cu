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

template <int scope, typename F>
__device__ __host__
bool check_simplex_3dclt(
    int current_timestep,
    const lattice4_t& domain, 
    const lattice4_t& core, 
    const lattice3_t& ext, // array dimension
    const element42_t& e,
    const F *uv[2],
    cp_t &p)
{
  if (e.corner[3] != current_timestep)
    return false;
  
  int vertices[3][4], indices[3];
  size_t local_indices[3];
  for (int i = 0; i < 3; i ++) {
    for (int j = 0; j < 4; j ++) {
      vertices[i][j] = e.corner[j]
        + unit_simplex_offset_4_2<scope>(e.type, i, j);
      if (vertices[i][j] < domain.st[j] || 
          vertices[i][j] > domain.st[j] + domain.sz[j] - 1)
        return false;
    }
    indices[i] = domain.to_index(vertices[i]);
    local_indices[i] = ext.to_index(vertices[i]);
  }
 
  F X[3][4], A[3][3];
  F u[3], v[3];
  for (int i = 0; i < 3; i ++) {
    const size_t k = local_indices[i]; // k = ext.to_index(vertices[i]);
    const size_t t = unit_simplex_offset_4_2<scope>(e.type, i, 3);
      
    u[i] = uv[t][k*2];
    v[i] = uv[t][k*2+1];

    for (int j = 0; j < 3; j ++)
      X[i][j] = vertices[i][j] * h[0]->cell_lengths[j] + h[0]->origins[j];
    X[i][3] = vertices[i][3];
  }

  // locate zero
  F mu[3], // barycentric coordinates
    cond; // condition number
  inverse_lerp_s2v2(psi, mu, &cond);

  // interpolation
  F x[4];
  lerp_s2v4(X, mu, x);

  // result
  p.x[0] = x[0];
  p.x[1] = x[1];
  p.x[2] = x[2];
  p.t = x[3];
  p.cond = cond;

  return true;
}

template <int scope, typename F>
__global__
void sweep_simplices(
    int current_timestep,
    const lattice4_t domain,
    const lattice4_t core,
    const lattice3_t ext, // array dimension
    const F *uv_c, // current timestep
    const F *uv_n, // next timestep
    unsigned long long &ncps, cp_t *cps)
{
  const F *uv[2] = {uv_c, uv_n};
  
  int tid = getGlobalIdx_3D_1D();
  const element42_t e = element42_from_index<scope>(core, tid);
  
  cp_t cp;
  bool succ = check_simplex_3dclt<scope>(
      current_timestep,
      domain, core, ext, e, uv, cp);

  if (succ) {
    unsigned long long i = atomicAdd(&ncps, 1ul);
    cp.tag = tid;
    cps[i] = cp;
  }
}

template <int scope>
static std::vector<cp_t> extract_3dclt(
    int current_timestep,
    const lattice4_t& domain,
    const lattice4_t& core, 
    const lattice3_t& ext,
    const double *uv_c,
    const double *uv_n)
{
  auto t0 = std::chrono::high_resolution_clock::now();

  const size_t ntasks = core.n() * ntypes_4_2<scope>();
  // fprintf(stderr, "ntasks=%zu\n", ntasks);
  const int maxGridDim = 1024;
  const int blockSize = 256;
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;
  
  if (nBlocks >= maxGridDim) 
    gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else 
    gridSize = dim3(nBlocks);

  F *duv_c = NULL, *duv_n = NULL;
  if (duv_c) {
    cudaMalloc((void**)&_c, sizeof(double) * ext.n());
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: allocating _c");
    cudaMemcpy(_c, uv_c, sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: copying _c");
  }
  if (duv_n) {
    cudaMalloc((void**)&_n, sizeof(double) * ext.n());
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: allocating _l");
    cudaMemcpy(_n, uv_n, sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: copying _l");
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
      domain, core, ext, duv_c, duv_n, 
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
 
  if (duv_c) cudaFree(duv_c);
  if (duv_n) cudaFree(duv_n);
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
extract_3dclt_cuda(
    int scope, 
    int current_timestep, 
    const ftk::lattice& domain,
    const ftk::lattice& core, 
    const ftk::lattice& ext, 
    const double *uv_c, 
    const double *uv_l)
{
  lattice4_t D(domain);
  lattice4_t C(core);
  lattice3_t E(ext);

  if (scope == scope_interval) 
    return extract_3dclt<scope_interval>(current_timestep, D, C, E, uv_c, uv_l);
  else
    return extract_3dclt<scope_ordinal>(current_timestep, D, C, E, uv_c, uv_l);
}
