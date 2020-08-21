#ifndef _FTK_LINEAR_INTERPOLATION_H
#define _FTK_LINEAR_INTERPOLATION_H

#include <ftk/ftk_config.hh>

namespace ftk {

template <typename T>
__device__ __host__
inline T lerp_s1(const T v[2], const T mu[2])
{
  return v[0] * mu[0] + v[1] * mu[1];
}

template <typename T>
__device__ __host__
inline void lerp_s1v2(const T V[2][2], const T mu[2], T v[2])
{
  v[0] = V[0][0] * mu[0] + V[1][0] * mu[1];
  v[1] = V[0][1] * mu[0] + V[1][1] * mu[1];
}

template <typename T>
__device__ __host__
inline void lerp_s1v2(const T V[2][2], const T mu, T v[2])
{
  const T nu = T(1) - mu;
  v[0] = V[0][0] * mu + V[1][0] * nu;
  v[1] = V[0][1] * mu + V[1][1] * nu;
}

template <typename T>
__device__ __host__
inline void lerp_s1v3(const T V[2][3], const T mu[2], T v[3])
{
  v[0] = V[0][0] * mu[0] + V[1][0] * mu[1];
  v[1] = V[0][1] * mu[0] + V[1][1] * mu[1];
  v[2] = V[0][2] * mu[0] + V[1][2] * mu[1];
}

template <typename T>
__device__ __host__
inline T lerp_s2(const T v[3], const T mu[3])
{
  return v[0] * mu[0] + v[1] * mu[1] + v[2] * mu[2];
}

template <typename T, int n>
inline void lerp_s2v(const T V[3][n], const T mu[3], T v[n])
{
  for (int i = 0; i < n; i ++)
    v[i] = V[0][i] * mu[0] + V[1][i] * mu[1] + V[2][i] * mu[2];
}

template <typename T>
inline void lerp_s2v2(const T V[3][2], const T mu[3], T v[2])
{
  v[0] = V[0][0] * mu[0] + V[1][0] * mu[1] + V[2][0] * mu[2];
  v[1] = V[0][1] * mu[0] + V[1][1] * mu[1] + V[2][1] * mu[2];
}

template <typename T>
__device__ __host__
inline void lerp_s2v3(const T V[3][3], const T mu[3], T v[3])
{
  v[0] = V[0][0] * mu[0] + V[1][0] * mu[1] + V[2][0] * mu[2];
  v[1] = V[0][1] * mu[0] + V[1][1] * mu[1] + V[2][1] * mu[2];
  v[2] = V[0][2] * mu[0] + V[1][2] * mu[1] + V[2][2] * mu[2];
}

template <typename T>
__device__ __host__
inline void lerp_s2m2x2(const T V[3][2][2], const T mu[3], T v[2][2])
{
  v[0][0] = V[0][0][0] * mu[0] + V[1][0][0] * mu[1] + V[2][0][0] * mu[2];
  v[0][1] = V[0][0][1] * mu[0] + V[1][0][1] * mu[1] + V[2][0][1] * mu[2];
  v[1][0] = V[0][1][0] * mu[0] + V[1][1][0] * mu[1] + V[2][1][0] * mu[2];
  v[1][1] = V[0][1][1] * mu[0] + V[1][1][1] * mu[1] + V[2][1][1] * mu[2];
}

template <typename T>
__device__ __host__
inline T lerp_s3(const T v[4], const T mu[4])
{
  return v[0] * mu[0] + v[1] * mu[1] + v[2] * mu[2] + v[3] * mu[3];
}

template <typename T>
__device__ __host__
inline void lerp_s3v3(const T V[4][3], const T mu[4], T v[3])
{
  v[0] = V[0][0] * mu[0] + V[1][0] * mu[1] + V[2][0] * mu[2] + V[3][0] * mu[3];
  v[1] = V[0][1] * mu[0] + V[1][1] * mu[1] + V[2][1] * mu[2] + V[3][1] * mu[3];
  v[2] = V[0][2] * mu[0] + V[1][2] * mu[1] + V[2][2] * mu[2] + V[3][2] * mu[3];
}

template <typename T>
__device__ __host__
inline void lerp_s3v4(const T V[4][4], const T mu[4], T v[4])
{
  v[0] = V[0][0] * mu[0] + V[1][0] * mu[1] + V[2][0] * mu[2] + V[3][0] * mu[3];
  v[1] = V[0][1] * mu[0] + V[1][1] * mu[1] + V[2][1] * mu[2] + V[3][1] * mu[3];
  v[2] = V[0][2] * mu[0] + V[1][2] * mu[1] + V[2][2] * mu[2] + V[3][2] * mu[3];
  v[3] = V[0][3] * mu[0] + V[1][3] * mu[1] + V[2][3] * mu[2] + V[3][3] * mu[3];
}

template <typename T>
__device__ __host__
inline void lerp_s3m3x3(const T V[4][3][3], const T mu[4], T v[3][3])
{
  for (int j = 0; j < 3; j ++)
    for (int k = 0; k < 3; k ++) {
      v[j][k] = T(0);
      for (int i = 0; i < 4; i ++) 
        v[j][k] += V[i][j][k] * mu[i];
    }
}

}

#endif
