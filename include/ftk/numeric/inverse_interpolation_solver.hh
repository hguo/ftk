#ifndef _FTK_INVERSE_INTERPOLATION_SOLVER_HH
#define _FTK_INVERSE_INTERPOLATION_SOLVER_HH

namespace ftk {

template <typename T>
bool inverse_linear_interpolation1(ValueType v[2], ValueType lambda[1])
{
  return 0;
}

template <typename T>
bool inverse_barycentric3_interpolation2(ValueType v[3][2], ValueType lambda[3])
{
  return false;
}

bool inverse_barycentric4_interpolation3(ValueType v[4][3], ValueType lambda[4])
{
  return false;
}

bool inverse_bilinear2_interpolation2(ValueType v[4][2], ValueType lambda[2])
{
  return false;
}

bool inserve_barycentric3_linear_interpolation3(ValueType v[3][3], ValueType, w[3][3], ValueType lambda[4])
{
  return false;
}

}

#endif
