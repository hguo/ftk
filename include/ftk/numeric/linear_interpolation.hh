#ifndef _FTK_LINEAR_INTERPOLATION_H
#define _FTK_LINEAR_INTERPOLATION_H

namespace ftk {

template <typename T>
inline T linear_interpolation2(const T v[2], const T mu[2])
{
  return v[0] * mu[0] + v[1] * mu[1];
}

template <typename T>
inline void linear_interpolation2_2(const T V[2][2], const T mu[2], T v[2])
{
  v[0] = V[0][0] * mu[0] + V[1][0] * mu[1];
  v[1] = V[0][1] * mu[0] + V[1][1] * mu[1];
}

template <typename T>
inline T linear_interpolation3(const T v[3], const T mu[3])
{
  return v[0] * mu[0] + v[1] * mu[1] + v[2] * mu[2];
}

template <typename T>
inline void linear_interpolation3_3(const T V[9], const T mu[3], T v[3])
{
  v[0] = V[0] * mu[0] + V[1] * mu[1] + V[2] * mu[2];
  v[1] = V[3] * mu[0] + V[4] * mu[1] + V[5] * mu[2];
  v[2] = V[6] * mu[0] + V[7] * mu[1] + V[8] * mu[2];
}

template <typename T>
inline void linear_interpolation3_3(const T V[3][3], const T mu[3], T v[3])
{
  v[0] = V[0][0] * mu[0] + V[1][0] * mu[1] + V[2][0] * mu[2];
  v[1] = V[0][1] * mu[0] + V[1][1] * mu[1] + V[2][1] * mu[2];
  v[2] = V[0][2] * mu[0] + V[1][2] * mu[1] + V[2][2] * mu[2];
}

template <typename T>
inline void linear_interpolation4_3(const T V[4][3], const T mu[4], T v[3])
{
  v[0] = V[0][0] * mu[0] + V[1][0] * mu[1] + V[2][0] * mu[2] + V[3][0] * mu[3];
  v[1] = V[0][1] * mu[0] + V[1][1] * mu[1] + V[2][1] * mu[2] + V[3][1] * mu[3];
  v[2] = V[0][2] * mu[0] + V[1][2] * mu[1] + V[2][2] * mu[2] + V[3][2] * mu[3];
}

template <typename T>
inline void linear_interpolation3_2x2(const T V[3][2][2], const T mu[3], T v[2][2])
{
  v[0][0] = V[0][0][0] * mu[0] + V[1][0][0] * mu[1] + V[2][0][0] * mu[2];
  v[0][1] = V[0][0][1] * mu[0] + V[1][0][1] * mu[1] + V[2][0][1] * mu[2];
  v[1][0] = V[0][1][0] * mu[0] + V[1][1][0] * mu[1] + V[2][1][0] * mu[2];
  v[1][1] = V[0][1][1] * mu[0] + V[1][1][1] * mu[1] + V[2][1][1] * mu[2];
}

}

#endif
