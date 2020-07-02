#ifndef _FTK_COND_HH
#define _FTK_COND_HH

#include <ftk/numeric/matrix_inverse.hh>
#include <ftk/numeric/matrix_norm.hh>

namespace ftk {

template <typename T>
T cond_real2x2(const T A[2][2])
{
  T invA[2][2];
  matrix_inverse2x2(A, invA);

  // return matrix_2norm_real2x2(invA) * matrix_2norm_real2x2(A);
  return matrix_1norm_real2x2<T>(invA) * matrix_1norm_real2x2(A);
}

template <typename T>
T cond_real3x3(const T A[3][3])
{
  T invA[3][3];
  matrix_inverse3x3(A, invA);

  // return matrix_2norm_real2x2(invA) * matrix_2norm_real2x2(A);
  return matrix_1norm_real3x3<T>(invA) * matrix_1norm_real3x3(A);
}

}

#endif
