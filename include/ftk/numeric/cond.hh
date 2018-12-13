#ifndef _FTK_COND_HH
#define _FTK_COND_HH

#include <ftk/numeric/matrix_inverse.hh>
#include <ftk/numeric/matrix_norm.hh>

namespace ftk {

template <typename T>
T cond_real2x2(const T A[2][2])
{
  T invA[2][2];
  matrix_inverse2(A, invA);

  // return matrix_2norm_real2x2(invA) * matrix_2norm_real2x2(A);
  return matrix_1norm_real2x2(invA) * matrix_1norm_real2x2(A);
}

}

#endif
