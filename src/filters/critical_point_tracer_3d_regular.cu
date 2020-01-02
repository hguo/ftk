#include <nvfunctional>
#include <cstdio>
#include <cassert>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/hypermesh/lattice.hh>
#include <ftk/filters/critical_point.hh>
#include "common.cuh"

__device__
bool detect_cp3t(
    const double v[4][3],
    const int vertices[4][4], 
    cp4_t &cp)
{
  double mu[4];
  bool succ = ftk::inverse_lerp_s3v3(v, mu);
 
  if (succ) {
    double X[4][4];
    for (int i = 0; i < 4; i ++)
      for (int j = 0; j < 4; j ++)
        X[i][j] = vertices[i][j];
    ftk::lerp_s3v4(X, mu, cp.x);
    return true;
  } else 
    return false;
}

template <int scope>
__device__
bool check_simplex_cp3t(
    int current_timestep,
    const lattice4_t& domain, 
    const lattice4_t& core, 
    const lattice3_t& ext, // array dimension
    const element43_t& e, 
    const double *V[2],
    cp4_t &cp)
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

  return detect_cp3t(v, vertices, cp);
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
    unsigned long long &ncps, cp4_t *cps)
{
  const double *V[2] = {Vl, Vc};
  
  int tid = getGlobalIdx_3D_1D();
  const element43_t e = element43_from_index<scope>(core, tid);

  cp4_t cp;
  bool succ = check_simplex_cp3t<scope>(
      current_timestep,
      domain, core, ext, e, V, cp);

  if (succ) {
    unsigned long long i = atomicAdd(&ncps, 1ul);
    cp.tag = tid;
    cps[i] = cp;
  }
}

template <int scope>
static std::vector<cp4_t> extract_cp3dt(
    int current_timestep,
    const lattice4_t& domain,
    const lattice4_t& core, 
    const lattice3_t& ext, 
    const double *Vc, 
    const double *Vl)
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

  unsigned long long *dncps; // number of cps
  cudaMalloc((void**)&dncps, sizeof(unsigned long long));
  cudaMemset(dncps, 0, sizeof(unsigned long long));

  cp4_t *dcps;
  cudaMalloc((void**)&dcps, sizeof(cp4_t) * core.n());
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMalloc/cudaMemcpy");

  fprintf(stderr, "calling kernel func...\n");
  sweep_simplices<scope><<<gridSize, blockSize>>>(
      current_timestep, 
      domain, core, ext, dVc, dVl, 
      *dncps, dcps);
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices, kernel function");

  unsigned long long ncps = 0;
  cudaMemcpy(&ncps, dncps, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMemcpyDeviceToHost, dncps");
  fprintf(stderr, "ncps=%lu\n", ncps);

  std::vector<cp4_t> cps(ncps);
  cudaMemcpy(cps.data(), dcps, sizeof(cp4_t) * ncps, cudaMemcpyDeviceToHost);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMemcpyDeviceToHost");
  
  cudaFree(dVc);
  cudaFree(dVl);
  cudaFree(dncps);
  cudaFree(dcps);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaFree");
 
  cudaDeviceSynchronize();
  fprintf(stderr, "exitting gpu kernel, ncps=%lu\n", ncps);

  return cps;
}

std::vector<cp4_t>
extract_cp3dt_cuda(
    int scope, 
    int current_timestep, 
    const ftk::lattice& domain,
    const ftk::lattice& core, 
    const ftk::lattice& ext, 
    const double *Vc, 
    const double *Vl)
{
  lattice4_t D(domain);
  lattice4_t C(core);
  lattice3_t E(ext);

  if (scope == scope_interval) 
    return extract_cp3dt<scope_interval>(current_timestep, D, C, E, Vc, Vl);
  if (scope == scope_ordinal) 
    return extract_cp3dt<scope_ordinal>(current_timestep, D, C, E, Vc, Vl);
  else // scope == 2
    return extract_cp3dt<scope_all>(current_timestep, D, C, E, Vc, Vl);
}
