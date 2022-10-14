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
    const int nc,
    const F *uv[2],
    cp_t &p)
{
  // typedef ftk::fixed_point<> fp_t;
  
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
 
  F X[3][4], UV[3][2];
  // fp_t UVf[3][2];
  const F factor = 32768;
  int64_t UVf[3][2];
  for (int i = 0; i < 3; i ++) {
    const size_t k = local_indices[i]; // k = ext.to_index(vertices[i]);
    const size_t t = unit_simplex_offset_4_2<scope>(e.type, i, 3);
    
    UV[i][0] = uv[t][k*nc];
    UV[i][1] = uv[t][k*nc+1];

    UVf[i][0] = UV[i][0] * factor;
    UVf[i][1] = UV[i][1] * factor;

    for (int j = 0; j < 3; j ++)
      X[i][j] = vertices[i][j]; 
    X[i][3] = vertices[i][3];
  }

  bool succ = robust_critical_point_in_simplex2(UVf, indices);
  if (succ) {
    // locate zero
    F mu[3], // barycentric coordinates
      cond; // condition number
    inverse_lerp_s2v2(UV, mu, &cond);
    ftk::clamp_barycentric<3>(mu);

    // interpolation
    F x[4];
    lerp_s2v4(X, mu, x);

    // result
    p.x[0] = x[0];
    p.x[1] = x[1];
    p.x[2] = x[2];
    p.t = x[3];
    // p.cond = cond;

    // printf("%f, %f, %f, %f\n", x[0], x[1], x[2], x[3]);

    return true;
  } else 
    return false;
}

template <int scope, typename F>
__global__
void sweep_simplices(
    int current_timestep,
    const lattice4_t domain,
    const lattice4_t core,
    const lattice3_t ext, // array dimension
    const int nc,
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
      domain, core, ext, e, nc, uv, cp);

  if (succ) {
    unsigned long long i = atomicAdd(&ncps, 1ul);
    cp.tag = tid;
    cps[i] = cp;
  }
}

template <int scope, typename F>
static std::vector<cp_t> extract_3dclt(
    int current_timestep,
    const lattice4_t& domain,
    const lattice4_t& core, 
    const lattice3_t& ext,
    const int nc,
    const float *uv_c,
    const float *uv_n)
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
  if (uv_c) {
    cudaMalloc((void**)&duv_c, sizeof(float) * ext.n() * nc);
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: allocating uv_c");
    cudaMemcpy(duv_c, uv_c, sizeof(float) * ext.n() * nc, cudaMemcpyHostToDevice);
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: copying uv_c");
  }
  if (uv_n) {
    cudaMalloc((void**)&duv_n, sizeof(float) * ext.n() * nc);
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: allocating uv_l");
    cudaMemcpy(duv_n, uv_n, sizeof(float) * ext.n() * nc, cudaMemcpyHostToDevice);
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: copying uv_l");
  }
  
  unsigned long long *dncps; // number of cps
  cudaMalloc((void**)&dncps, sizeof(unsigned long long));
  cudaMemset(dncps, 0, sizeof(unsigned long long));
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: allocating dncps");

  cp_t *dcps;
  cudaMalloc((void**)&dcps, 1024*1024*256); // sizeof(cp_t) * core.n());
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: allocating dcps");
  cudaDeviceSynchronize();

  fprintf(stderr, "calling kernel func...\n");
  sweep_simplices<scope><<<gridSize, blockSize>>>(
      current_timestep, 
      domain, core, ext, nc, duv_c, duv_n, 
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
  float duration = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() * 1e-9;

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
    const int nc,
    const float *uv_c, 
    const float *uv_n)
{
  lattice4_t D(domain);
  lattice4_t C(core);
  lattice3_t E(ext);

  if (scope == scope_interval) 
    return extract_3dclt<scope_interval, float>(current_timestep, D, C, E, nc, uv_c, uv_n);
  else
    return extract_3dclt<scope_ordinal, float>(current_timestep, D, C, E, nc, uv_c, uv_n);
}
