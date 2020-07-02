#ifndef _UTILS_CUH
#define _UTILS_CUH

#include <cstdio>

inline void checkCuda(cudaError_t e, const char *situation) {
  if (e != cudaSuccess) {
    fprintf(stderr, "CUDA Error: %s: %s\n", situation, cudaGetErrorString(e));
  }
}

inline void checkLastCudaError(const char *situation) {
  // cudaDeviceSynchronize();
  checkCuda(cudaGetLastError(), situation);
}

static inline int idivup(int a, int b)
{
  return (a%b!=0) ? (a/b+1) : (a/b); 
}

static inline dim3 idivup(dim3 a, dim3 b)
{
  return dim3(idivup(a.x, b.x), idivup(a.y, b.y), idivup(a.z, b.z));
}

template <typename T>
__device__
inline static T fmod1(T x, T y)
{
  T z = fmod(x, y);
  if (z<0) z += y;
  return z;
}

template <typename T>
__device__
inline static T mod2pi(T x)
{
  T y = fmod(x, 2*M_PI); 
  if (y<0) y+= 2*M_PI;
  return y; 
}

template <typename T>
__device__
inline static T mod2pi1(T x)
{
  return mod2pi(x + M_PI) - M_PI;
}

template <typename T> 
__device__
inline int sgn(T x) 
{
  return (T(0) < x) - (x < T(0));
}

template <typename T>
__device__
static inline T inner_product(const T A[3], const T B[3])
{
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

template <typename T>
__device__
static inline T dist2(const T A[3], const T B[3])
{
  const T D[3] = {B[0]-A[0], B[1]-A[1], B[2]-A[2]};
  return inner_product(D, D);
}


#endif
