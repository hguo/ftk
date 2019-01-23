#ifndef _FTK_VORTEX_CRITERIA_HH
#define _FTK_VORTEX_CRITERIA_HH

#include <ftk/numeric/matrix_decomposition.hh>
#include <ftk/numeric/matrix_multiplication.hh>
#include <ftk/numeric/matrix_addition.hh>
#include <ftk/numeric/eigen_solver.hh>
#include <complex>
#include <algorithm>

namespace ftk {

template <class T>
inline T vortex_lambda2_criterion(const T J[3][3])
{
  T S[3][3], O[3][3];
  matrix_symmetric_decomposition_3x3(J, S);
  matrix_antisymmetric_decomposition_3x3(J, O);

  T squareS[3][3], squareO[3][3];
  matrix_square_3x3(S, squareS);
  matrix_square_3x3(O, squareO);

  T M[3][3];
  matrix_addition3x3(squareS, squareO, M);

  std::complex<T> eig[3];
  std::complex<T> eigvec[3][3];
  eig3(M, eig, eigvec);

  std::sort(eig, eig+3, [](const std::complex<T>& a, const std::complex<T>& b) -> bool {
      return std::abs(a) > std::abs(b);
  });

  return std::abs(eig[2]);
}

}

#endif
