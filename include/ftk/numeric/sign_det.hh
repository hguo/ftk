#ifndef _FTK_SIGN_DET_HH
#define _FTK_SIGN_DET_HH

#include <ftk/config.hh>
#include <ftk/int128.hh>
#include <ftk/numeric/sign.hh>
#include <ftk/numeric/det.hh>
#include <ftk/numeric/swap.hh>

// reference:
// Edelsbrunner and Mucke, Simulation of simplicity: A technique to cope with degenerate cases in geometric algorithms.

namespace ftk {

template <typename I, typename T>
__device__ __host__
inline bool robust_smaller(I i, I j, I k, I l, T pij, T pkl)
{
  if (pij != pkl) return pij < pkl;
  else if (i != k) return i > k;
  else return j < l;
}

template <typename T=int128_t>
__device__ __host__
inline int robust_sign_det2(const T X[2])
{
  for (int t = 0; t < 2; t ++) {
    int sigma = 0;
    if (t == 0) {
      const T M[2][2] = {
        {X[0], T(1)}, 
        {X[1], T(1)}
      };
      sigma = sign(det2(M));
    } else 
      sigma = 1;

    if (sigma != 0) return sigma;
  }
  // assert(false);
  return 0;
}

template <typename T=int128_t>
__device__ __host__
inline int robust_sign_det3(const T X[3][2])
{
  for (int t = 0; t < 5; t ++) {
    int sigma = 0;
    if (t == 0) {
      const T M[3][3] = {
        {X[0][0], X[0][1], T(1)}, 
        {X[1][0], X[1][1], T(1)},
        {X[2][0], X[2][1], T(1)}
      };
      // print3x3("M", M);
      // std::cerr << "det=" << det3(M) << std::endl;
      sigma = sign(det3(M));
    } else if (t == 1) {
      const T M[2][2] = {
        {X[1][0], T(1)},
        {X[2][0], T(1)}
      };
      // print2x2("M", M);
      // std::cerr << "det=" << det2(M) << std::endl;
      sigma = -sign(det2(M));
    } else if (t == 2) {
      const T M[2][2] = {
        {X[1][1], T(1)},
        {X[2][1], T(1)}
      };
      // print2x2("M", M);
      // std::cerr << "det=" << det2(M) << std::endl;
      sigma = sign(det2(M));
    } else if (t == 3) {
      const T M[2][2] = {
        {X[0][0], T(1)},
        {X[2][0], T(1)}
      };
      sigma = sign(det2(M));
    } else 
      sigma = 1;
      
    // fprintf(stderr, "t=%d, sigma=%d\n", t, sigma);
    if (sigma != 0)
      return sigma;
  }
  // assert(false);
  return 0; // useless
}

template <typename T=int128_t>
__device__ __host__
inline int robust_sign_det4(const T X[4][3])
{
  for (int t = 0; t < 15; t ++) {
    int sigma = 0;
    if (t == 0) {
      const T M[4][4] = {
        {X[0][0], X[0][1], X[0][2], T(1)},
        {X[1][0], X[1][1], X[1][2], T(1)},
        {X[2][0], X[2][1], X[2][2], T(1)},
        {X[3][0], X[3][1], X[3][2], T(1)}
      };
      // print4x4("M", M);
      // std::cerr << "det=" << det4(M) << std::endl;
      sigma = sign(det4(M));
    } else if (t == 1) {
      const T M[3][3] = {
        {X[1][0], X[1][1], T(1)},
        {X[2][0], X[2][1], T(1)},
        {X[3][0], X[3][1], T(1)}
      };
      sigma = sign(det3(M));
    } else if (t == 2) {
      const T M[3][3] = {
        {X[1][0], X[1][2], T(1)},
        {X[2][0], X[2][2], T(1)},
        {X[3][0], X[3][2], T(1)}
      };
      sigma = -sign(det3(M));
    } else if (t == 3) {
      const T M[3][3] = {
        {X[1][1], X[1][2], T(1)},
        {X[2][1], X[2][2], T(1)},
        {X[3][1], X[3][2], T(1)}
      };
      sigma = sign(det3(M));
    } else if (t == 4) {
      const T M[3][3] = {
        {X[0][0], X[0][1], T(1)},
        {X[2][0], X[2][1], T(1)},
        {X[3][0], X[3][1], T(1)}
      };
      sigma = -sign(det3(M));
    } else if (t == 5) {
      const T M[2][2] = {
        {X[2][0], T(1)},
        {X[3][0], T(1)}
      };
      sigma = sign(det2(M));
    } else if (t == 6) {
      const T M[2][2] = {
        {X[2][1], T(1)},
        {X[3][1], T(1)}
      };
      sigma = -sign(det2(M));
    } else if (t == 7) {
      const T M[3][3] = {
        {X[0][0], X[0][2], T(1)},
        {X[2][0], X[2][2], T(1)},
        {X[3][0], X[3][2], T(1)}
      };
      sigma = sign(det3(M));
    } else if (t == 8) {
      const T M[2][2] = {
        {X[2][2], T(1)},
        {X[3][2], T(1)}
      };
      sigma = sign(det2(M));
    } else if (t == 9) {
      const T M[3][3] = {
        {X[0][1], X[0][2], T(1)},
        {X[2][1], X[2][2], T(1)},
        {X[3][1], X[3][2], T(1)}
      };
      sigma = -sign(det3(M));
    } else if (t == 10) {
      const T M[3][3] = {
        {X[0][0], X[0][1], T(1)},
        {X[1][0], X[1][1], T(1)},
        {X[3][0], X[3][1], T(1)}
      };
      sigma = sign(det3(M));
    } else if (t == 11) {
      const T M[2][2] = {
        {X[1][0], T(1)},
        {X[3][0], T(1)}
      };
      sigma = -sign(det2(M));
    } else if (t == 12) {
      const T M[2][2] = {
        {X[1][1], T(1)},
        {X[3][1], T(1)}
      };
      sigma = sign(det2(M));
    } else if (t == 13) {
      const T M[2][2] = {
        {X[0][0], T(1)},
        {X[3][0], T(1)}
      };
      sigma = sign(det2(M));
    } else 
      sigma = 1;

    if (sigma != 0) return sigma;
  }
  // assert(false);
  return 0; // useless
}

// returns number of swaps for bubble sort
template <int n, typename T=int128_t>
__device__ __host__
inline int nswaps_bubble_sort(T arr[n], T order[n])
{
  for (int i = 0; i < n; i ++)
    order[i] = i;

  int nswaps = 0;
  for (int i = 0; i < n-1; i ++) 
    for (int j = 0; j < n-i-1; j ++) 
      if (arr[j] > arr[j+1]) {
        swap_helper(arr[j], arr[j+1]);
        swap_helper(order[j], order[j+1]);
        nswaps ++;
      }
  
  return nswaps;
}

template <typename T=int128_t, typename I=int128_t>
__device__ __host__
inline I positive1(const T X1[2], const I indices1[2])
{
  I indices[2], orders[2];
  for (int i = 0; i < 2; i ++)
    indices[i] = indices1[i];
  int s = nswaps_bubble_sort<2, I>(indices, orders);

  T X[2];
  for (int i = 0; i < 2; i ++)
    X[i] = X1[orders[i]];

  I d = robust_sign_det2(X);

  if (s % 2 != 0) 
    d = -d;

  return d;
}

template <typename T=int128_t, typename I=int128_t>
__device__ __host__
inline I positive2(const T X1[3][2], const I indices1[3])
{
  I indices[3], orders[3];
  for (int i = 0; i < 3; i ++)
    indices[i] = indices1[i];
  int s = nswaps_bubble_sort<3, I>(indices, orders); // number of swaps to get sorted indices
  // fprintf(stderr, "nswaps=%d\n", s);

  T X[3][2];
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 2; j ++)
      X[i][j] = X1[orders[i]][j];
      // X[i][j] = X1[i][j];

  // fprintf(stderr, "indices=%d, %d, %d\n", indices[0], indices[1], indices[2]);
  // print3x2("X", X);
  int d = robust_sign_det3(X);

  if (s % 2 != 0) // odd 
    d = -d;

  return d;
}

template <typename T=int128_t, typename I=int128_t>
__device__ __host__
inline I positive3(const T X1[4][3], const I indices1[4])
{
  I indices[4], orders[4];
  for (int i = 0; i < 4; i ++)
    indices[i] = indices1[i];
  I s = nswaps_bubble_sort<4, I>(indices, orders);

  T X[4][3];
  for (int i = 0; i < 4; i ++)
    for (int j = 0; j < 3; j ++)
      X[i][j] = X1[orders[i]][j];

  I d = robust_sign_det4(X);

  if (s % 2 != 0) // odd
    d = -d;

  return d;
}

template <typename T>
__device__ __host__
inline bool robust_intersect_half_line2(
    const T v[2], 
    int j, const T vj[2], 
    int k, const T vk[2])
{
  // assuming 
  if (!robust_smaller(j, 1, k, 1, vj[1], vk[1]))
    return robust_intersect_half_line2(v, k, vk, j, vj);
  
  if (robust_smaller(j, 1, 0, 1, vj[1], v[1]) && 
      robust_smaller(0, 1, k, 1, v[1], vk[1]))
  {
    // assuming 0 < i < j
    const T M[3][2] = {
      {v[0], v[1]},
      {vj[0], vj[1]},
      {vk[0], vk[1]}
    };
    const int d = robust_sign_det3(M);
    fprintf(stderr, "d=%d\n", d);
    return d == 1; 
  } else 
    return false;
}

template <typename T>
__device__ __host__
inline bool robust_point_in_polygon2(
    const T x[2], const int n, const int indices[], const T v[][2])
{
  int count = 0;
  for (int i = 0; i < n-1; i ++) {
    if (robust_intersect_half_line2(x, indices[i], v[i], indices[i+1], v[i+1]))
      count ++;
  }
  fprintf(stderr, "count=%d\n", count);
  if (count % 2 == 1) return true;
  else return false;
}

// check if a point in a 1-simplex
template <typename T=long long>
__device__ __host__
inline bool robust_point_in_simplex1(const T X[2], const int indices[2], const T x, const int ix)
{
  const int s = positive1(X, indices);

  for (int i = 0; i < 2; i ++) {
    T Y[2];
    int my_indices[2];
    for (int j = 0; j < 2; j ++)
      if (i == j) {
        my_indices[j] = ix;
        Y[j] = x;
      } else {
        my_indices[j] = indices[j];
        Y[j] = X[j];
      }

    int si = positive1(Y, my_indices);
    if (s != si) 
      return false;
  }
  return true;
}

// check if a point is in a 2-simplex
template <typename T=int128_t, typename I=int128_t>
__device__ __host__
inline bool robust_point_in_simplex2(const T X[3][2], const I indices[3], const T x[2], const I ix) //, const int sign=1)
{
  // print3x2("X", X);
  const I s = positive2(X, indices); // orientation of the simplex
  // fprintf(stderr, "orientation s=%d\n", s);
  for (int i = 0; i < 3; i ++) {
    T Y[3][2];
    I my_indices[3];
    for (int j = 0; j < 3; j ++)
      if (i == j) {
        my_indices[j] = ix;
        for (int k = 0; k < 2; k ++) 
          Y[j][k] = x[k];
      } else {
        my_indices[j] = indices[j];
        for (int k = 0; k < 2; k ++) 
          Y[j][k] = X[j][k];
      }
  
    // print3x2("Y", Y);

    I si = positive2(Y, my_indices);
    // fprintf(stderr, "s=%d, s[%d]=%d\n", s, i, si);
    if (s != si)
      return false;
  }
  return true;
}

template <typename T=int128_t, typename I=int128_t>
__device__ __host__
inline bool robust_point_in_simplex3(const T X[4][3], const I indices[3], const T x[3], I ix)
{
  I s = positive3(X, indices);
  for (int i = 0; i < 4; i ++) {
    T Y[4][3];
    I my_indices[4];
    for (int j = 0; j < 4; j ++)
      if (i == j) {
        my_indices[j] = ix;
        for (int k = 0; k < 3; k ++) 
          Y[j][k] = x[k];
      } else {
        my_indices[j] = indices[j];
        for (int k = 0; k < 3; k ++) 
          Y[j][k] = X[j][k];
      }

    I si = positive3(Y, my_indices);
    if (s != si) return false;
  }
  return true;
}

} // namespace

#endif
