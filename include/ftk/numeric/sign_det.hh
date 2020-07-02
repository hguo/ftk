#ifndef _FTK_SIGN_DET_HH
#define _FTK_SIGN_DET_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/det.hh>

// reference:
// Edelsbrunner and Mucke, Simulation of simplicity: A technique to cope with degenerate cases in geometric algorithms.

namespace ftk {

template <typename T>
__device__ __host__
inline bool robust_smaller(int i, int j, int k, int l, T pij, T pkl)
{
  if (pij != pkl) return pij < pkl;
  else if (i != k) return i > k;
  else return j < l;
}

template <typename T=long long>
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

template <typename T=long long>
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

template <typename T=long long>
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

template <typename T>
__device__ __host__
inline void swap_helper(T& a, T& b)
{
  T c(a); a = b; b = c;
}

// returns number of swaps for bubble sort
template <int n, typename T>
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

template <typename T=long long>
__device__ __host__
inline int positive2(const T X1[3][2], const int indices1[3])
{
  int indices[3], orders[3];
  for (int i = 0; i < 3; i ++)
    indices[i] = indices1[i];
  int s = nswaps_bubble_sort<3, int>(indices, orders); // number of swaps to get sorted indices
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

template <typename T=long long>
__device__ __host__
inline int positive3(const T X1[4][3], const int indices1[4])
{
  int indices[4], orders[4];
  for (int i = 0; i < 3; i ++)
    indices[i] = indices1[i];
  int s = nswaps_bubble_sort<4, int>(indices, orders);

  T X[4][3];
  for (int i = 0; i < 4; i ++)
    for (int j = 0; j < 3; j ++)
      X[i][j] = X1[orders[i]][j];

  int d = robust_sign_det4(X);

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
  // assume j < k
  if (k > j)
    return robust_intersect_half_line2(v, k, vk, j, vj);

  if (robust_smaller(j, 1, 0, 1, vj[1], v[1]) && 
      robust_smaller(0, 1, k, 1, v[1], vk[1]))
  {
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

// check if a point is in a 2-simplex
template <typename T=long long>
__device__ __host__
inline bool robust_point_in_simplex2(const T X[3][2], const int indices[3], const T x[2], const int ix) //, const int sign=1)
{
  // print3x2("X", X);
  const int s = positive2(X, indices); // orientation of the simplex
  // fprintf(stderr, "orientation s=%d\n", s);
  for (int i = 0; i < 3; i ++) {
    T Y[3][2];
    int my_indices[3];
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

    int si = positive2(Y, my_indices);
    // fprintf(stderr, "s=%d, s[%d]=%d\n", s, i, si);
    if (s != si)
      return false;
  }
  return true;
}

template <typename T=long long>
__device__ __host__
inline bool robust_point_in_simplex3(const T X[4][3], const int indices[3], const T x[3], int ix)
{
  int s = positive3(X, indices);
  for (int i = 0; i < 4; i ++) {
    T Y[4][3];
    for (int j = 0; j < 4; j ++)
      if (i == j)
        for (int k = 0; k < 3; k ++) 
          Y[j][k] = x[k];
      else 
        for (int k = 0; k < 3; k ++) 
          Y[j][k] = X[j][k];

    int si = positive3(Y, indices);
    if (s != si) return false;
  }
  return true;
}

} // namespace

#endif
