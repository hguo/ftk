#ifndef _FTK_EIGEN_SOLVER2_HH
#define _FTK_EIGEN_SOLVER2_HH

#include <ftk/numeric/trace.hh>
#include <ftk/numeric/det.hh>
#include <ftk/numeric/matrix_multiplication.hh>
#include <ftk/numeric/quadratic_solver.hh>
#include <ftk/numeric/cubic_solver.hh>
#include <ftk/numeric/quartic_solver.hh>
#include <ftk/numeric/vector_normalization.hh>
#include <ftk/numeric/characteristic_polynomial.hh>
#include <ftk/numeric/print.hh>
#include <iostream>

namespace ftk {

template <typename T>
__device__ __host__
inline void solve_eigenvalues_symmetric2x2(T m00, T m10, T m11, T eig[2])
{
  const T b = -(m00 + m11), c = m00*m11 - m10*m10;
  // const T delta = b*b - 4*c;
#ifdef __CUDACC__
  const T delta = __fma_rn(b, b, -4*c); // TODO: __fma_rn is only for doubles
  const T sqrt_delta = delta < 0 ? 0 : sqrt(delta); // in theory, delta should not be less than 0
#else
  const T delta = std::fma(b, b, -4*c); 
  const T sqrt_delta = delta < 0 ? 0 : std::sqrt(delta); // in theory, delta should not be less than 0
#endif
  
  eig[0] = T(0.5) * (-b+sqrt_delta);
  eig[1] = T(0.5) * (-b-sqrt_delta);

  if (abs(eig[0]) < abs(eig[1])) {
#ifdef __CUDACC__
    T t(eig[0]);
    eig[0] = eig[1];
    eig[1] = t;
#else
    std::swap(eig[0], eig[1]);
#endif
  }
}

template <typename T>
inline void solve_eigenvalues_symmetric2x2(const T m[2][2], T eig[2])
{
  return solve_eigenvalues_symmetric2x2(m[0][0], m[1][0], m[1][1], eig);
}

template <typename T>
__device__ __host__
inline int solve_eigenvalues2x2(const T M[2][2], T eig[2])
{
  T P[3];
  characteristic_polynomial_2x2(M, P);
  // print2x2("M", M);
  // fprintf(stderr, "P=%f, %f, %f\n", P[0], P[1], P[2]);
  return solve_quadratic_real(P, eig);
}

template <typename T>
__device__ __host__
inline T solve_eigenvalues2x2(const T M[2][2], std::complex<T> eig[2])
{
  T P[3];
  characteristic_polynomial_2x2(M, P);
  return solve_quadratic(P, eig); // returns delta
}

template <typename T>
inline void solve_eigenvectors2x2(const T M[2][2], int n, const T eig[2], T eigvecs[2][2])
{
  for (int i = 0; i < n; i ++) {
    const T a[2] = {
      M[0][0] - eig[i] - M[0][1], 
      M[1][0] - M[1][1] + eig[i]
    };
    const T b[2] = {-M[0][1], -M[1][1] + eig[i]};

    solve_least_square2x1(a, b, eigvecs[i][0]);
    eigvecs[i][1] = T(1) - eigvecs[i][0];
    vector_normalization2_2(eigvecs[i]);
  }
}

template <typename T>
inline int solve_generalized_eigenvalues2x2(const T A[2][2], const T B[2][2], T eig[2])
{
  T P[3];
  characteristic_polynomial_2x2(A, B, P);
  return solve_quadratic_real(P, eig);
}

}
#endif
