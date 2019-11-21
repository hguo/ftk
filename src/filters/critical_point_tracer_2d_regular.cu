#include <nvfunctional>
#include <cstdio>
// #include <ftk/filters/critical_point_tracker_2d.hh>
#include "threadIdx.cuh"
#include "utils.cuh"
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/hypermesh/lattice.hh>
#include <ftk/filters/critical_point.hh>

template <int N=3>
struct lite_lattice_t {
  int st[N], sz[N], prod[N];

  __host__ lite_lattice_t(const ftk::lattice& L) {
    for (int i = 0; i < N; i ++) {
      st[i] = L.starts_[i];
      sz[i] = L.sizes_[i];
      prod[i] = L.prod_[i];
    }
  };

  lite_lattice_t(int st_[N], int sz_[N]) {
    for (int i = 0; i < N; i ++) {
      st[i] = st_[i];
      sz[i] = sz_[i];
      if (i == 0) prod[i] = 1;
      else prod[i] = prod[i-1] * sz[i-1];
    }
  }

  template <typename uint=size_t>
  __device__ __host__
  uint n() const {return prod[N-1] * sz[N-1];}

  template <typename uint=size_t>
  __device__ __host__
  uint to_index(const int idx1[N]) const {
    int idx[N];
    for (int j = 0; j < N; j ++)
      idx[j] = idx1[j] - st[j];

    uint i = idx[0];
    for (int j = 1; j < N; j ++)
      i += idx[j] * prod[j];

    return i;
  }

  template <typename uint=size_t>
  __device__ __host__
  void from_index(uint i, int idx[N]) const {
    for (int j = N-1; j > 0; j --) {
      idx[j] = i / prod[j];
      i -= idx[j] * prod[j];
    }
    idx[0] = i;

    for (int j = 0; j < N; j ++)
      idx[j] += st[j];
  }
};

template <typename T, int N=3>
T access(const lite_lattice_t<N>& l, const T *p, const int idx[N])
{
  size_t i = l.to_index(idx);
  return p[i];
}

template <int N=3>
struct lite_element_t {
  int corner[N], /*d,*/ type;
};

typedef lite_lattice_t<3> lattice3_t;
typedef lite_element_t<3> element32_t;
typedef ftk::critical_point_t<3, double> cp3_t;
  
__device__ __constant__ 
int ntypes_3[4] = {1, 7, 12, 6}, 
    // unit_simplices_3_3[6][4] = {0}, // TODO
    unit_simplices_3_2[12][3][3] = {
      {{0,0,0},{0,0,1},{0,1,1}},
      {{0,0,0},{0,0,1},{1,0,1}},
      {{0,0,0},{0,0,1},{1,1,1}},
      {{0,0,0},{0,1,0},{0,1,1}},
      {{0,0,0},{0,1,0},{1,1,0}},
      {{0,0,0},{0,1,0},{1,1,1}},
      {{0,0,0},{0,1,1},{1,1,1}},
      {{0,0,0},{1,0,0},{1,0,1}},
      {{0,0,0},{1,0,0},{1,1,0}},
      {{0,0,0},{1,0,0},{1,1,1}},
      {{0,0,0},{1,0,1},{1,1,1}},
      {{0,0,0},{1,1,0},{1,1,1}}
    };
    // unit_simplices_3_1[7][2] = {0};
  
template <typename uint=size_t>
__device__ __host__
element32_t element32_from_index(const lattice3_t& l, int scope, uint i) {
  element32_t e; // TODO
  
  e.type = i % 12; // m.ntypes(dim, scope);
  uint ii = i / 12; // m.ntypes(dim, scope);
  l.from_index(ii, e.corner);

  return e;
}

template <typename uint=size_t>
__device__ __host__
uint element32_to_index(const lattice3_t& l, int scope, const int idx[3]) {
  size_t i = l.to_index(idx);
  return i * 12; // m.ntypes(dim, scope);
}
  
__device__
bool check_simplex_cp2t(
    const lattice3_t& core, 
    const lattice3_t& ext, 
    const element32_t& e, 
    const double *V, 
    cp3_t &cp)
{
  int vertices[3][3];
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 3; j ++) {
      vertices[i][j] = e.corner[j]
        + unit_simplices_3_2[e.type][i][j];
      if (vertices[i][j] < core.st[j] || 
          vertices[i][j] > core.st[j] + core.sz[j] - 1)
        return false;
    }

  double v[3][2];
  for (int i = 0; i < 3; i ++) {
    size_t k = ext.to_index(vertices[i]);
    for (int j = 0; j < 2; j ++)
      v[i][j] = V[k*2+j]; // V has two channels
  }

  double mu[3];
  bool succ = ftk::inverse_lerp_s2v2(v, mu, 0.0);
 
  if (succ) {
    double X[3][3];
    for (int i = 0; i < 3; i ++)
      for (int j = 0; j < 3; j ++)
        X[i][j] = vertices[i][j];
    ftk::lerp_s2v3(X, mu, cp.x);
    return true;
  } else 
    return false;
}

__global__
void sweep_simplices(
    const lattice3_t core, int scope,
    const lattice3_t ext, const double *V, 
    unsigned long long &ncps, cp3_t *cps)
{
  int tid = getGlobalIdx_3D_1D();
  const element32_t e = element32_from_index(core, scope, tid);

  cp3_t cp;
  bool succ = check_simplex_cp2t(core, ext, e, V, cp);
  if (succ) {
    unsigned long long i = atomicAdd(&ncps, 1ul);
    cp.tag = tid;
    cps[i] = cp;
  }
}

static void extract_cp2dt(
    const lattice3_t& core, int scope, 
    const lattice3_t& ext, const double *V/* 4D array: 2*W*H*T */)
{
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
  cudaMalloc((void**)&dV, 2 * sizeof(double) * ext.n());
  cudaMemcpy(dV, V, 2 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);

  unsigned long long *dncps; // number of cps
  cudaMalloc((void**)&dncps, sizeof(unsigned long long));
  cudaMemset(dncps, 0, sizeof(unsigned long long));

  cp3_t *dcps;
  cudaMalloc((void**)&dcps, sizeof(cp3_t) * ext.n());
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMalloc/cudaMemcpy");

  fprintf(stderr, "calling kernel func...\n");
  sweep_simplices<<<gridSize, blockSize>>>(core, scope, ext, dV, *dncps, dcps);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices");

  unsigned long long ncps;
  cudaMemcpy(&ncps, dncps, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
  fprintf(stderr, "ncps=%lu\n", ncps);

  std::vector<cp3_t> cps(ncps);
  cudaMemcpy(cps.data(), dcps, sizeof(cp3_t) * ncps, cudaMemcpyDeviceToHost);
  
  cudaFree(dV);
  cudaFree(dncps);
  cudaFree(dcps);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaFree");
 
  cudaDeviceSynchronize();
  fprintf(stderr, "exit, ncps=%lu\n", ncps);
}

void extract_cp2dt_cuda(
    const ftk::lattice& core, int scope, 
    const ftk::lattice& ext, const double *V)
{
  lattice3_t C(core), E(ext);
  extract_cp2dt(C, scope, E, V);
}
