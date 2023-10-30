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
inline T solve_least_square4x3(const T A[4][3], const T b[4], T x[3], 
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

  matrix3x4_vector4_multiplication(invATAAT, b, x);
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

//////
template <typename T>
__device__ __host__
void elgs(int n, T **A, int *indx)
{
  for (int i = 0; i < n; i ++)
    indx[i] = i;

  // rescaling factors
  T c[n];
  for (int i = 0; i < n; i ++) {
    T c1 = 0.0;
    for (int j = 0; j < n; j ++)
      c1 = std::max(c1, std::abs(A[i][j]));
    c[i] = c1;
  }

  // search for the pivoting element
  T pi, pj;
  int k;
  for (int j = 0; j < n-1; j ++) {
    T pi1 = 0.0;
    for (int i = j; i < n; i ++) {
      pi = std::abs(A[indx[i]][j]) / c[indx[i]];
      if (pi > pi1) {
        pi1 = pi;
        k = i;
      }
    }
    
    // interchange rows
    std::swap(indx[j], indx[k]);
    for (int i = j; i < n; i ++) {
      pj = A[indx[i]][j] / A[indx[j]][j];
      A[indx[i]][j] = pj;

      for (int k = j; k < n; k ++) 
        A[indx[i]][k] = A[indx[i]][k] - pj * A[indx[j]][k];
    }
  }
}

template <typename T>
__device__ __host__
void legs(int n, T **A, T *b, T *x, int *indx)
{
  elgs(n, A, indx);

  for (int i = 0; i < n-1; i ++)
    for (int j = i; j < n; j ++)
      b[indx[j]] = b[indx[j]] - A[indx[j]][i] * b[indx[i]];

  x[n-1] = b[indx[n-1]] / A[indx[n-1]][n-1];
  for (int i = n-2; i >= 0; i --) {
    x[i] = b[indx[i]];
    for (int j = i; j < n; j ++)
      x[i] = x[i] - A[indx[i]][j] * x[j];
    x[i] = x[i] / A[indx[i]][i];
  }
}

}

#endif
