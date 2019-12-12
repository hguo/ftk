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

typedef lite_lattice_t<3> lattice3_t;
typedef lite_element_t<3> element32_t;
typedef ftk::critical_point_t<3, double> cp3_t;

typedef lite_lattice_t<4> lattice4_t;
typedef lite_element_t<4> element43_t;
typedef ftk::critical_point_t<4, double> cp4_t;
  
__device__ __constant__ 
int ntypes_3[4] = {1, 7, 12, 6}, 
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

__device__ __constant__
int ordinal_unit_simplex_indices_4_3[] = {
  16, 20, 30, 34, 46, 50
};

__device__ __constant__
int interval_unit_simplex_indices_4_3[] = {
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 
  13, 14, 15, 17, 18, 19, 21, 22, 23, 24, 25, 
  26, 27, 28, 29, 31, 32, 33, 35, 36, 37, 38, 
  39, 40, 41, 42, 43, 44, 45, 47, 48, 49, 51, 
  52, 53, 54, 55, 56, 57, 58, 59
};

__device__ __constant__ 
int unit_simplices_4_3[][4][4] = {
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

template <int scope>
__device__ __host__ inline int ntypes_3_2();

template <>
__device__ __host__ inline int ntypes_3_2<0>() { return 12; }

template <>
__device__ __host__ inline int ntypes_3_2<1>() { return 2; }

template <>
__device__ __host__ inline int ntypes_3_2<2>() { return 10; }

template <int scope>
__device__ __host__ inline int ntypes_4_3();

template <>
__device__ __host__ inline int ntypes_4_3<0>() { return 60; }

template <>
__device__ __host__ inline int ntypes_4_3<1>() { return 6; }

template <>
__device__ __host__ inline int ntypes_4_3<2>() { return 54; }

template <int scope=0, typename uint=size_t>
__device__ __host__
element32_t element32_from_index(const lattice3_t& l, uint i) {
  element32_t e; // TODO
  
  e.type = i % ntypes_3_2<scope>();
  uint ii = i / ntypes_3_2<scope>();
  l.from_index(ii, e.corner);

  return e;
}

template <int scope=0, typename uint=size_t>
__device__ __host__
uint element32_to_index(const lattice3_t& l, const int idx[3]) {
  size_t i = l.to_index(idx);
  return i * ntypes_3_2<scope>();
}

template <int scope=0, typename uint=size_t>
__device__ __host__
element43_t element43_from_index(const lattice4_t& l, uint i) {
  element43_t e;
  
  e.type = i % ntypes_4_3<scope>();
  uint ii = i / ntypes_4_3<scope>(); 
  l.from_index(ii, e.corner);

  return e;
}

template <int scope=0, typename uint=size_t>
__device__ __host__
uint element43_to_index(const lattice4_t& l, const int idx[4]) {
  size_t i = l.to_index(idx);
  return i * ntypes_4_3<scope>; 
}

#endif
