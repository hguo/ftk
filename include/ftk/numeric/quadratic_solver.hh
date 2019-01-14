#ifndef _FTK_QUADRATIC_SOLVER_HH
#define _FTK_QUADRATIC_SOLVER_HH

#include <cmath>
#include <complex>
#include <ftk/numeric/sqrt.hh>
#include <ftk/numeric/linear_solver.hh>

namespace ftk {

template <typename T>
inline void solve_quadratic(const T P[3], std::complex<T> x[2])
{
  const T delta = P[1]*P[1] - 4*P[2]*P[0];
  if (delta >= 0) {
    x[0] = (-P[1] + sqrt(delta)) / (2 * P[2]);
    x[1] = (-P[1] - sqrt(delta)) / (2 * P[2]);
  } else {
    x[0] = (-P[1] + complex_sqrt<T>(delta)) / (2 * P[2]);
    x[1] = (-P[1] - complex_sqrt<T>(delta)) / (2 * P[2]);
  }
}

template <typename T>
inline int solve_quadratic_real(const T P[3], T x[2])
{
  if (P[2] == T(0) || std::isinf(P[2]) || std::isnan(P[2])) return solve_linear_real1(P, x);

  const T delta = P[1]*P[1] - 4*P[2]*P[0];
  if (delta > 0) {
    x[0] = (-P[1] + sqrt(delta)) / (2 * P[2]);
    x[1] = (-P[1] - sqrt(delta)) / (2 * P[2]);
    return 2;
  } else if (delta == T(0)) {
    x[0] = -P[1] / (2 * P[2]);
    return 1;
  } else 
    return 0;
}

}

#endif
