#ifndef _FTK_LINEAR_SOLVE_H
#define _FTK_LINEAR_SOLVE_H

#include <ftk/numeric/det.hh>
#include <ftk/numeric/cond.hh>
#include <ftk/numeric/inner_product.hh>
#include <ftk/numeric/matrix_inverse.hh>
#include <ftk/numeric/print.hh>

namespace ftk {

template <typename T>
__device__ __host__
inline T solve_linear3x3(const T A[3][3], const T b[3], T x[3])
{
  T invA[3][3];
  const T det = matrix_inverse3x3(A, invA);
  matrix3x3_vector3_multiplication(invA, b, x);
  return det;
}

template <typename T> // returns determinant
__device__ __host__
inline T solve_linear2x2(const T A[2][2], const T b[2], T x[2], const T epsilon = std::numeric_limits<T>::epsilon())
{
  const T D  = det2(A),
          Dx = b[0]   *A[1][1] - A[0][1]*b[1],
          Dy = A[0][0]*b[1]    - b[0]   *A[1][0];

  x[0] = Dx / D; 
  x[1] = Dy / D;
  
  return D;
}

template <typename T>
__device__ __host__
inline T solve_least_square3x2(const T A[3][2], const T b[3], T x[2], const T epsilon = std::numeric_limits<T>::epsilon())
{
  T AT[2][3];
  transpose3x2(A, AT);

  T ATA[2][2];
  matrix2x3_matrix3x2_multiplication(AT, A, ATA);

  T invATA[2][2];
  const T det = matrix_inverse2x2(ATA, invATA);
  const T cond = cond_real2x2(ATA);

  T invATAAT[2][3];
  matrix2x2_matrix2x3_multiplication(invATA, AT, invATAAT);

  matrix2x3_vector3_multiplication(invATAAT, b, x);
  return cond;
}

template <typename T>
__device__ __host__
inline T solve_least_square3x2_2(const T A[3][2], const T B[3][2], T x[2][2], 
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  // print3x2("A", A);
  T AT[2][3];
  transpose3x2(A, AT);

  T ATA[2][2];
  matrix2x3_matrix3x2_multiplication(AT, A, ATA);
  // print2x2("ATA", ATA);

  T invATA[2][2];
  const T det = matrix_inverse2x2(ATA, invATA);
  const T cond = cond_real2x2(ATA);
  // print2x2("invATA", invATA);

  T invATAAT[2][3];
  matrix2x2_matrix2x3_multiplication(invATA, AT, invATAAT);

  matrix2x3_matrix3x2_multiplication(invATAAT, B, x);
  return cond;
}

template <typename T>
__device__ __host__
inline T solve_least_square4x3_3(const T A[4][3], const T B[4][3], T x[3][3], 
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  // print3x2("A", A);
  T AT[3][4];
  transpose4x3(A, AT);

  T ATA[3][3];
  matrix3x4_matrix4x3_multiplication(AT, A, ATA);
  // print2x2("ATA", ATA);

  T invATA[3][3];
  const T det = matrix_inverse3x3(ATA, invATA);
  const T cond = cond_real3x3(ATA);
  // print2x2("invATA", invATA);

  T invATAAT[3][4];
  matrix3x3_matrix3x4_multiplication(invATA, AT, invATAAT);

  matrix3x4_matrix4x3_multiplication(invATAAT, B, x);
  return cond;
}

}

#endif
