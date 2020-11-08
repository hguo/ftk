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
bool check_simplex_tdgl_vortex_3dt(
    int current_timestep,
    const lattice4_t& domain, 
    const lattice4_t& core, 
    const lattice3_t& ext, // array dimension
    const element43_t& e, 
    const double *Rho[2], // current and next timesteps
    const double *Phi[2], 
    cp_t &cp)
{
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
  
  double rho[3], phi[3];
  for (int i = 0; i < 3; i ++) {
    const size_t k = local_indices[i]; // k = ext.to_index(vertices[i]);
    const size_t t = unit_simplex_offset_4_3<scope>(e.type, i, 3);
      
    rho[i] = Rho[t][k];
    phi[i] = Phi[t][k];
  }

  return false; // WIP
}

template <int scope>
__global__
void sweep_simplices(
    int current_timestep,
    const lattice4_t domain,
    const lattice4_t core,
    const lattice3_t ext, // array dimension
    const double *rho_c, // current timestep
    const double *rho_n, // next timestep
    const double *phi_c, 
    const double *phi_n,
    unsigned long long &ncps, cp_t *cps)
{
  const double *Rho[2] = {rho_c, rho_n};
  const double *Phi[2] = {phi_c, phi_n};
  
  int tid = getGlobalIdx_3D_1D();
  const element43_t e = element43_from_index<scope>(core, tid);

  cp_t cp;
  bool succ = check_simplex_tdgl_vortex_3dt<scope>(
      current_timestep,
      domain, core, ext, e, V, J, S, cp);

  if (succ) {
    unsigned long long i = atomicAdd(&ncps, 1ul);
    cp.tag = tid;
    cps[i] = cp;
  }
}

template <int scope>
static std::vector<cp_t> extract_tdgl_vortex_3dt(
    int current_timestep,
    const lattice4_t& domain,
    const lattice4_t& core, 
    const lattice3_t& ext, 
    const double *rho_c,
    const double *rho_n, 
    const double *phi_c,
    const double *phi_n)
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

  double *drho_c = NULL, *drho_n = NULL;
  if (rho_c) {
    cudaMalloc((void**)&drho_c, 3 * sizeof(double) * ext.n());
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: allocating drho_c");
    cudaMemcpy(drho_c, rho_c, 3 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: copying drho_c");
  }
  if (rho_n) {
    cudaMalloc((void**)&drho_n, 3 * sizeof(double) * ext.n());
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: allocating drho_l");
    cudaMemcpy(drho_n, rho_n, 3 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
    checkLastCudaError("[FTK-CUDA] error: sweep_simplices: copying drho_l");
  }
  
  double *dphi_c = NULL, *dphi_n = NULL;
  if (phi_c) {
    cudaMalloc((void**)&dphi_c, 9 * sizeof(double) * ext.n());
    cudaMemcpy(dphi_c, phi_c, 9 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }
  if (phi_n) {
    cudaMalloc((void**)&dphi_n, 9 * sizeof(double) * ext.n());
    cudaMemcpy(dphi_n, phi_n, 9 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
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
      domain, core, ext, drho_c, drho_n, dphi_c, dphi_n, 
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
  
  if (drho_c) cudaFree(drho_c);
  if (drho_n) cudaFree(drho_n);
  if (dphi_c) cudaFree(dphi_c);
  if (dphi_n) cudaFree(dphi_n);
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
extract_tdgl_vortex_3dt_cuda(
    int scope, 
    int current_timestep, 
    const ftk::lattice& domain,
    const ftk::lattice& core, 
    const ftk::lattice& ext, 
    const double *rho_c, 
    const double *rho_l,
    const double *phi_c, 
    const double *phi_l)
{
  lattice4_t D(domain);
  lattice4_t C(core);
  lattice3_t E(ext);

  if (scope == scope_interval) 
    return extract_tdgl_vortex_3dt<scope_interval>(current_timestep, D, C, E, rho_c, rho_l, phi_c, phi_l);
  if (scope == scope_ordinal) 
    return extract_tdgl_vortex_3dt<scope_ordinal>(current_timestep, D, C, E, rho_c, rho_l, phi_c, phi_l);
  else // scope == 2
    return extract_tdgl_vortex_3dt<scope_all>(current_timestep, D, C, E, rho_c, rho_l, phi_c, phi_l);
}