#include <nvfunctional>
#include <cstdio>
// #include <ftk/filters/critical_point_tracker_2d.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/hypermesh/lattice.hh>
#include <ftk/filters/critical_point.hh>
#include "common.cuh"

typedef lite_lattice_t<4> lattice4_t;
typedef lite_element_t<4> element43_t;
typedef ftk::critical_point_t<4, double> cp4_t;
  
__device__ __constant__ 
int unit_simplices_4_3[60][4][4] = {
  {{0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 1}, {0, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 1}, {1, 0, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 0, 1}, {0, 1, 0, 1}, {0, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 0, 1}, {0, 1, 0, 1}, {1, 1, 0, 1}},
  {{0, 0, 0, 0}, {0, 0, 0, 1}, {0, 1, 0, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 0, 1}, {0, 1, 1, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 0, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 1, 0, 1}},
  {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, 1, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 1, 0, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 1, 1}, {0, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 1, 1}, {1, 0, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 1, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 1, 0}, {0, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 1, 0}, {1, 1, 1, 0}},
  {{0, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 1, 0}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 1, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 1, 0}, {1, 0, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 1, 0}, {1, 1, 1, 0}},
  {{0, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 1, 0}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 1, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 1, 0}, {1, 1, 1, 0}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 1, 1}, {0, 1, 1, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 0, 1, 1}, {1, 0, 1, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 1, 0, 0}, {0, 1, 0, 1}, {0, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 1, 0, 0}, {0, 1, 0, 1}, {1, 1, 0, 1}},
  {{0, 0, 0, 0}, {0, 1, 0, 0}, {0, 1, 0, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 1, 0, 0}, {0, 1, 1, 0}, {0, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 1, 0, 0}, {0, 1, 1, 0}, {1, 1, 1, 0}},
  {{0, 0, 0, 0}, {0, 1, 0, 0}, {0, 1, 1, 0}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 1, 0, 0}, {0, 1, 1, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 1, 0, 0}, {1, 1, 0, 0}, {1, 1, 0, 1}},
  {{0, 0, 0, 0}, {0, 1, 0, 0}, {1, 1, 0, 0}, {1, 1, 1, 0}},
  {{0, 0, 0, 0}, {0, 1, 0, 0}, {1, 1, 0, 0}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 1, 0, 0}, {1, 1, 0, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 1, 0, 0}, {1, 1, 1, 0}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 1, 0, 1}, {0, 1, 1, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 1, 0, 1}, {1, 1, 0, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 1, 1, 0}, {0, 1, 1, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {0, 1, 1, 0}, {1, 1, 1, 0}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 0, 1}, {1, 0, 1, 1}},
  {{0, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 0, 1}, {1, 1, 0, 1}},
  {{0, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 0, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 1, 0}, {1, 0, 1, 1}},
  {{0, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 1, 0}, {1, 1, 1, 0}},
  {{0, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 1, 0}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 1, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 0}, {1, 1, 0, 1}},
  {{0, 0, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 0}, {1, 1, 1, 0}},
  {{0, 0, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 0}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {1, 0, 0, 0}, {1, 1, 1, 0}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {1, 0, 0, 1}, {1, 0, 1, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {1, 0, 0, 1}, {1, 1, 0, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {1, 0, 1, 0}, {1, 0, 1, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {1, 0, 1, 0}, {1, 1, 1, 0}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {1, 1, 0, 0}, {1, 1, 0, 1}, {1, 1, 1, 1}},
  {{0, 0, 0, 0}, {1, 1, 0, 0}, {1, 1, 1, 0}, {1, 1, 1, 1}}
};
  
template <typename uint=size_t>
__device__ __host__
element43_t element43_from_index(const lattice4_t& l, int scope, uint i) {
  element43_t e; // TODO
  
  e.type = i % 60; // m.ntypes(dim, scope);
  uint ii = i / 60; // m.ntypes(dim, scope);
  l.from_index(ii, e.corner);

  return e;
}

template <typename uint=size_t>
__device__ __host__
uint element43_to_index(const lattice4_t& l, int scope, const int idx[4]) {
  size_t i = l.to_index(idx);
  return i * 60; // m.ntypes(dim, scope);
}
  
__device__
bool check_simplex_cp3t(
    const lattice4_t& core, 
    const lattice4_t& ext, 
    const element43_t& e, 
    const double *V, 
    cp4_t &cp)
{
  int vertices[4][4];
  for (int i = 0; i < 4; i ++)
    for (int j = 0; j < 4; j ++) {
      vertices[i][j] = e.corner[j]
        + unit_simplices_4_3[e.type][i][j];
      if (vertices[i][j] < core.st[j] || 
          vertices[i][j] > core.st[j] + core.sz[j] - 1)
        return false;
    }

  double v[4][3];
  for (int i = 0; i < 4; i ++) {
    size_t k = ext.to_index(vertices[i]);
    for (int j = 0; j < 3; j ++)
      v[i][j] = V[k*3+j]; // V has three channels
  }

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

__global__
void sweep_simplices(
    const lattice4_t core, int scope,
    const lattice4_t ext, const double *V, 
    unsigned long long &ncps, cp4_t *cps)
{
  int tid = getGlobalIdx_3D_1D();
  const element43_t e = element43_from_index(core, scope, tid);

  cp4_t cp;
  bool succ = check_simplex_cp3t(core, ext, e, V, cp);
  if (succ) {
    unsigned long long i = atomicAdd(&ncps, 1ul);
    cp.tag = tid;
    cps[i] = cp;
  }
}

static std::vector<cp4_t> extract_cp3dt(
    const lattice4_t& core, int scope, 
    const lattice4_t& ext, const double *V/* 5D array: 2*W*H*D*T */)
{
  fprintf(stderr, "init GPU...\n");
  const size_t ntasks = core.n() * 12; // ntypes_3[2] = 12; ntypes_3 is in device constant memory
  const int maxGridDim = 1024;
  const int blockSize = 256;
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;
  
  if (nBlocks >= maxGridDim) 
    gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else 
    gridSize = dim3(nBlocks);

  double *dV;
  cudaMalloc((void**)&dV, 3 * sizeof(double) * ext.n());
  cudaMemcpy(dV, V, 3 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);

  unsigned long long *dncps; // number of cps
  cudaMalloc((void**)&dncps, sizeof(unsigned long long));
  cudaMemset(dncps, 0, sizeof(unsigned long long));

  cp4_t *dcps;
  cudaMalloc((void**)&dcps, sizeof(cp4_t) * ext.n());
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMalloc/cudaMemcpy");

  fprintf(stderr, "calling kernel func...\n");
  sweep_simplices<<<gridSize, blockSize>>>(core, scope, ext, dV, *dncps, dcps);
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices, kernel function");

  unsigned long long ncps = 0;
  cudaMemcpy(&ncps, dncps, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMemcpyDeviceToHost, dncps");
  fprintf(stderr, "ncps=%lu\n", ncps);

  std::vector<cp4_t> cps(ncps);
  cudaMemcpy(cps.data(), dcps, sizeof(cp4_t) * ncps, cudaMemcpyDeviceToHost);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMemcpyDeviceToHost");
  
  cudaFree(dV);
  cudaFree(dncps);
  cudaFree(dcps);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaFree");
 
  cudaDeviceSynchronize();
  fprintf(stderr, "exit, ncps=%lu\n", ncps);

  return cps;
}

std::vector<cp4_t>
extract_cp3dt_cuda(
    const ftk::lattice& core, int scope, 
    const ftk::lattice& ext, const double *V)
{
  lattice4_t C(core), E(ext);
  return extract_cp3dt(C, scope, E, V);
}
