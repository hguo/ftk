#ifndef _FTK_LINEAR_INTERPOLATION_H
#define _FTK_LINEAR_INTERPOLATION_H

namespace ftk {

template <typename ValueType>
inline void linear_interpolation3(const ValueType V[9], const ValueType lambda[3], ValueType v[3])
{
  v[0] = V[0] * lambda[0] + V[1] * lambda[1] + V[2] * lambda[2];
  v[1] = V[3] * lambda[0] + V[4] * lambda[1] + V[5] * lambda[2];
  v[2] = V[6] * lambda[0] + V[7] * lambda[1] + V[8] * lambda[2];
}

template <typename ValueType>
inline void linear_interpolation3(const ValueType V[3][3], const ValueType lambda[3], ValueType v[3])
{
  v[0] = V[0][0] * lambda[0] + V[1][0] * lambda[1] + V[2][0] * lambda[2];
  v[1] = V[0][1] * lambda[0] + V[1][1] * lambda[1] + V[2][1] * lambda[2];
  v[2] = V[0][2] * lambda[0] + V[1][2] * lambda[1] + V[2][2] * lambda[2];
}

}

#endif
