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

template <typename ValueType>
inline void bilinear_interpolation3_location(const ValueType X[12], ValueType alpha, ValueType beta, ValueType x[3])
{
  ValueType p[2] = {alpha, beta};
  ValueType a[3], b[3]; 

  a[0] = (1-p[0])*X[0] + p[0]*X[1];
  a[1] = (1-p[0])*X[4] + p[0]*X[5];
  a[2] = (1-p[0])*X[8] + p[0]*X[9];

  b[0] = (1-p[0])*X[3] + p[0]*X[2];
  b[1] = (1-p[0])*X[7] + p[0]*X[6];
  b[2] = (1-p[0])*X[11]+ p[0]*X[10];

  x[0] = (1-p[1])*a[0] + p[1]*b[0];
  x[1] = (1-p[1])*a[1] + p[1]*b[1];
  x[2] = (1-p[1])*a[2] + p[1]*b[2];
}

}

#endif
