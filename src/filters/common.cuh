#ifndef _FTK_COMMON_CUH
#define _FTK_COMMON_CUH

#include "threadIdx.cuh"
#include "utils.cuh"

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

#endif
