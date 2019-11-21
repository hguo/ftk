#include <nvfunctional>
#include <cstdio>
// #include <ftk/filters/critical_point_tracker_2d.hh>
#include "threadIdx.cuh"
#include "utils.cuh"
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/hypermesh/lattice.hh>

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

template <int N=3>
struct lite_cp_t {
  size_t idx;
  double x[N]; 
  int type;
};
  
typedef lite_lattice_t<3> lattice3_t;
typedef lite_element_t<3> element32_t;
typedef lite_cp_t<3> cp3_t;
  
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
    int tid,
    const lattice3_t& domain, 
    const lattice3_t& block, 
    const element32_t& e, 
    double *V, 
    cp3_t &cp)
{
  int vertices[3][3];
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 3; j ++) {
      vertices[i][j] = e.corner[j]
        + unit_simplices_3_2[e.type][i][j];
      if (vertices[i][j] < block.st[j] || 
          vertices[i][j] > block.st[j] + block.sz[j] - 1)
        return false;
    }

  double v[3][2];
  for (int i = 0; i < 3; i ++) {
    size_t k = block.to_index(vertices[i]);
    for (int j = 0; j < 2; j ++)
      v[i][j] = V[k*2+j];
  }

  double mu[3];
  bool succ = ftk::inverse_lerp_s2v2(v, mu, 0.0);
 
  if (succ) {
    double X[3][3];
    for (int i = 0; i < 3; i ++)
      for (int j = 0; j < 3; j ++)
        X[i][j] = vertices[i][j];
    ftk::lerp_s2v3(X, mu, cp.x);
#if 0
    printf("tid=%d, mu=%f, %f, %f, v={%f, %f, %f, %f, %f, %f}, verts=(%d, %d, %d), (%d, %d, %d), (%d, %d, %d)\n",
        tid, mu[0], mu[1], mu[2], 
        v[0][0], v[0][1], v[1][0], v[1][1], v[2][0], v[2][1],
        vertices[0][0], vertices[0][1], vertices[0][2],
        vertices[1][0], vertices[1][1], vertices[1][2],
        vertices[2][0], vertices[2][1], vertices[2][2]);
#endif
    return true;
  } else 
    return false;
}

__global__
void sweep_simplices(
    const lattice3_t domain, int scope, 
    const lattice3_t block, double *V)
{
  int tid = getGlobalIdx_3D_1D();
  const element32_t e = element32_from_index(domain, scope, tid);

  cp3_t cp;
  bool succ = check_simplex_cp2t(tid, domain, block, e, V, cp);
#if 1
  if (succ && tid < 4000000)
    printf("succ, tid=%d, x=%f, %f, %f\n", tid, 
        cp.x[0], cp.x[1], cp.x[2]);
#endif
}

static void extract_cp2dt(
    const lattice3_t& domain, int scope, 
    const lattice3_t& block, double *V/* 4D array: 2*W*H*T */)
{
  const size_t ntasks = block.n() * 12; // ntypes_3[2] = 12; ntypes_3 is in device constant memory
  const int maxGridDim = 1024;
  const int blockSize = 256;
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;
  
  if (nBlocks >= maxGridDim) 
    gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else 
    gridSize = dim3(nBlocks);

  double *dV;
  cudaMalloc((void**)&dV, 2 * sizeof(double) * block.n());
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMalloc");

  cudaMemcpy(dV, V, 2*sizeof(double)*block.n(), cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMemcpy");

  fprintf(stderr, "calling kernel func...\n");
  sweep_simplices<<<gridSize, blockSize>>>(domain, scope, block, dV);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices");

  cudaFree(dV);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaFree");
  
  cudaDeviceSynchronize();
  fprintf(stderr, "exit.\n");
}

void extract_cp2dt(
    const ftk::lattice& domain, int scope, 
    const ftk::lattice& block, double *V)
{
  lattice3_t D(domain), B(block);
  std::cerr << domain << std::endl;
  std::cerr << block << std::endl;
  extract_cp2dt(D, scope, B, V);
}

#if 0
int main(int argc, char **argv)
{
  int st[3] = {0, 0, 0}, 
      sz[3] = {256, 256, 16};
  lattice3_t domain(st, sz);
  double *V = NULL;
  int scope = 0;

  extract_cp2dt(domain, scope, domain, V);
  return 0;
}
#endif




#if 0
struct cp2d_tracker_context {
public:
  __device__ __host__
  void element_for_3_2(
      const lattice3_t& l, int scope,
      size_t tid, // thread id
      const nvstd::function<bool(const element32_t&)> &f)
      // const nvstd::function<void()> &f)
  {
    element32_t e = element32_from_index(l, scope, tid);
    f(e);
    // f(e);
  }

};
#endif
