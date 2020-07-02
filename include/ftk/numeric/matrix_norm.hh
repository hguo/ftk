#ifndef _FTK_MATRIX_NORM_HH
#define _FTK_MATRIX_NORM_HH

#include <ftk/numeric/matrix_multiplication.hh>
#include <ftk/numeric/eigen_solver2.hh>
#include <ftk/numeric/eigen_solver3.hh>
#include <ftk/numeric/transpose.hh>

namespace ftk {

template <typename T>
inline T matrix_1norm_real2x2(const T m[2][2])
{
  // simply the maximum absolute column sum of the matrix (wikipedia)
  return std::max(std::abs(m[0][0]) + std::abs(m[1][0]), std::abs(m[0][1]) + std::abs(m[1][1]));
}

template <typename T>
inline T matrix_1norm_real3x3(const T m[3][3])
{
  // simply the maximum absolute column sum of the matrix (wikipedia)
  return std::max(
      std::max(
        std::abs(m[0][0]) + std::abs(m[1][0]) + std::abs(m[2][0]), 
        std::abs(m[0][1]) + std::abs(m[1][1]) + std::abs(m[2][1])), 
      std::abs(m[0][2]) + std::abs(m[1][2]) + std::abs(m[2][2]));
}

template <typename T>
inline T matrix_2norm_real2x2(const T A[2][2])
{
  T B[2][2]; // B is the transpose of A
  transpose2(A, B);
  
  T M[2][2]; // M = A^T * A
  matrix_matrix_multiplication_2x2_2x2(B, A, M);

  T eig[2];
  solve_eigenvalues_real_symmetric2(M, eig);

  return sqrt(eig[0]);
}

template <typename T>
inline T matrix_inf_norm_real2x2(const T m[2][2])
{
  // simply the maximum absolute row sum of the matrix (wikipedia)
  return std::max(std::abs(m[0][0]) + std::abs(m[0][1]), std::abs(m[1][0]) + std::abs(m[1][1]));
}

template <typename T>
inline T matrix_frobenius_norm_real2x2(const T m[2][2])
{
  return m[0][0]*m[0][0] + m[0][1]*m[0][1] + m[1][0]*m[1][0] + m[1][1]*m[1][1];
}

template <typename T>
inline T matrix_frobenius_norm_real2x3(const T m[2][3])
{
  return m[0][0]*m[0][0] + m[0][1]*m[0][1] + m[0][2]*m[0][2]
    + m[1][0]*m[1][0] + m[1][1]*m[1][1] + m[1][2]*m[1][2];
}

template <typename T>
inline T matrix_frobenius_norm_real3x2(const T m[3][2])
{
  return m[0][0]*m[0][0] + m[0][1]*m[0][1] 
    + m[1][0]*m[1][0] + m[1][1]*m[1][1]
    + m[2][0]*m[2][0] + m[2][1]*m[2][1];
}

}

#endif
