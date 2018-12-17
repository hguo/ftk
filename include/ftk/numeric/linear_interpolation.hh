#ifndef _FTK_LINEAR_INTERPOLATION_H
#define _FTK_LINEAR_INTERPOLATION_H

namespace ftk {

template <typename T>
inline void linear_interpolation3_3(const T V[9], const T lambda[3], T v[3])
{
  v[0] = V[0] * lambda[0] + V[1] * lambda[1] + V[2] * lambda[2];
  v[1] = V[3] * lambda[0] + V[4] * lambda[1] + V[5] * lambda[2];
  v[2] = V[6] * lambda[0] + V[7] * lambda[1] + V[8] * lambda[2];
}

template <typename T>
inline void linear_interpolation3_3(const T V[3][3], const T lambda[3], T v[3])
{
  v[0] = V[0][0] * lambda[0] + V[1][0] * lambda[1] + V[2][0] * lambda[2];
  v[1] = V[0][1] * lambda[0] + V[1][1] * lambda[1] + V[2][1] * lambda[2];
  v[2] = V[0][2] * lambda[0] + V[1][2] * lambda[1] + V[2][2] * lambda[2];
}

template <typename T>
inline void linear_interpolation4_3(const T V[4][3], const T lambda[4], T v[3])
{
  v[0] = V[0][0] * lambda[0] + V[1][0] * lambda[1] + V[2][0] * lambda[2] + V[3][0] * lambda[3];
  v[1] = V[0][1] * lambda[0] + V[1][1] * lambda[1] + V[2][1] * lambda[2] + V[3][1] * lambda[3];
  v[2] = V[0][2] * lambda[0] + V[1][2] * lambda[1] + V[2][2] * lambda[2] + V[3][2] * lambda[3];
}

template <typename T>
inline void linear_interpolation3_2x2(const T V[3][2][2], const T lambda[3], T v[2][2])
{
  v[0][0] = V[0][0][0] * lambda[0] + V[1][0][0] * lambda[1] + V[2][0][0] * lambda[2];
  v[0][1] = V[0][0][1] * lambda[0] + V[1][0][1] * lambda[1] + V[2][0][1] * lambda[2];
  v[1][0] = V[0][1][0] * lambda[0] + V[1][1][0] * lambda[1] + V[2][1][0] * lambda[2];
  v[1][1] = V[0][1][1] * lambda[0] + V[1][1][1] * lambda[1] + V[2][1][1] * lambda[2];
}

}

#endif
