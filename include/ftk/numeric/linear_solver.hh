#ifndef _FTK_LINEAR_SOLVE3_H
#define _FTK_LINEAR_SOLVE3_H

#include <ftk/numeric/det.hh>
#include <ftk/numeric/cond.hh>
#include <ftk/numeric/inner_product.hh>
#include <ftk/numeric/matrix_inverse.hh>
#include <ftk/numeric/print.hh>

namespace ftk {

template <typename T>
inline T linear_solver2_cond(const T A[2][2], const T b[2], T x[2])
{
  const T det = linear_solver2(A, b, x);
  const T cond = cond_real2x2(A);
  return cond;
}

template <typename T>
inline T linear_solver3(const T A[3][3], const T b[3], T x[3])
{
  T invA[3][3];
  const T det = matrix_inverse3x3(A, invA);
  matrix_vector_multiplication_3x3(invA, b, x);
  return det;
}

template <typename T> // returns determinant
// inline T linear_solver2(const T A[2][2], const T b[2], T x[2])
inline T solve_linear_real2x2(const T A[2][2], const T b[2], T x[2], const T epsilon = std::numeric_limits<T>::epsilon())
{
  const T D  = A[0][0]*A[1][1] - A[0][1]*A[1][0],
          Dx = b[0]   *A[1][1] - A[0][1]*b[1],
          Dy = A[0][0]*b[1]    - b[0]   *A[1][0];

  x[0] = Dx / D; 
  x[1] = Dy / D;
  // x[0] = (Dx + epsilon) / (D + epsilon);
  // x[1] = (Dy + epsilon) / (D + epsilon);
  
  return D;
}

template <typename T>
inline T solve_least_square3x2(const T A[3][2], const T b[3], T x[2], const T epsilon = std::numeric_limits<T>::epsilon())
{
  T AT[2][3];
  transpose3x2(A, AT);

  T ATA[2][2];
  matrix_matrix_multiplication_2x3_3x2(AT, A, ATA);

  T invATA[2][2];
  const T det = matrix_inverse2x2(ATA, invATA);
  const T cond = cond_real2x2(ATA);

  T invATAAT[2][3];
  matrix_matrix_multiplication_2x2_2x3(invATA, AT, invATAAT);

  matrix_vector_multiplication_2x3(invATAAT, b, x);
  return cond;
}

template <typename T>
inline T solve_least_square3x2_2(const T A[3][2], const T B[3][2], T x[2][2], const T epsilon = std::numeric_limits<T>::epsilon())
{
  // print3x2("A", A);
  T AT[2][3];
  transpose3x2(A, AT);

  T ATA[2][2];
  matrix_matrix_multiplication_2x3_3x2(AT, A, ATA);
  // print2x2("ATA", ATA);

  T invATA[2][2];
  const T det = matrix_inverse2x2(ATA, invATA);
  const T cond = cond_real2x2(ATA);
  // print2x2("invATA", invATA);

  T invATAAT[2][3];
  matrix_matrix_multiplication_2x2_2x3(invATA, AT, invATAAT);

  matrix_matrix_multiplication_2x3_3x2(invATAAT, B, x);
  return cond;
}

}

#endif
