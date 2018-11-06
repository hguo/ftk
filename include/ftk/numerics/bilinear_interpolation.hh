#ifndef _FTK_BILINEAR_INTERPOLATION_H
#define _FTK_BILINEAR_INTERPOLATION_H

namespace ftk {

template <typename ValueType>
inline void bilinear_interpolation3(const ValueType V[12], ValueType alpha, ValueType beta, ValueType v[3])
{
  const ValueType alpha1 = ValueType(1) - alpha, beta1 = ValueType(1) - beta;
  v[0] = V[0] * alpha1 * beta1 + V[1] * alpha * beta1 + V[3]  * alpha1 * beta + V[2]  * alpha * beta;
  v[1] = V[4] * alpha1 * beta1 + V[5] * alpha * beta1 + V[7]  * alpha1 * beta + V[6]  * alpha * beta;
  v[2] = V[8] * alpha1 * beta1 + V[9] * alpha * beta1 + V[11] * alpha1 * beta + V[10] * alpha * beta;
}

}

#endif
