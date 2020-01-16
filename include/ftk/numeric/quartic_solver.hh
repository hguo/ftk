#ifndef _FTK_QUARTIC_SOLVE_H
#define _FTK_QUARTIC_SOLVE_H

#include <cmath>
#include <complex>
#include <ftk/numeric/sqrt.hh>
#include <ftk/numeric/cbrt.hh>

namespace ftk {

// the followings are based on the code from https://github.com/sidneycadot/quartic/blob/master/solve-quartic.cc
// a * x^4 + b * x^3 + c * x^2 + d * x + e == 0
template <typename T>
void quartic_solve(T b_, T c_, T d_, T e_, std::complex<T> roots[4])
{
  // The algorithm below was derived by solving the quartic in Mathematica, and simplifying the resulting expression by hand.
  const std::complex<T> b(b_), c(c_), d(d_), e(e_);

  const std::complex<T> Q1 = c * c - T(3) * b * d + T(12) * e;
  const std::complex<T> Q2 = T(2) * c * c * c - T(9) * b * c * d + T(27) * d * d + T(27) * b * b * e - T(72) * c * e;
  const std::complex<T> Q3 = T(8) * b * c - T(16) * d - T(2) * b * b * b;
  const std::complex<T> Q4 = T(3) * b * b - T(8) * c;

  const std::complex<T> Q5 = complex_cbrt(Q2 / T(2) + complex_sqrt(Q2 * Q2 / T(4) - Q1 * Q1 * Q1));
  const std::complex<T> Q6 = (Q1 / Q5 + Q5) / T(3);
  const std::complex<T> Q7 = T(2) * complex_sqrt(Q4 / T(12) + Q6);

  roots[0] = (-b - Q7 - complex_sqrt(T(4) * Q4 / T(6) - T(4) * Q6 - Q3 / Q7)) / T(4);
  roots[1] = (-b - Q7 + complex_sqrt(T(4) * Q4 / T(6) - T(4) * Q6 - Q3 / Q7)) / T(4);
  roots[2] = (-b + Q7 - complex_sqrt(T(4) * Q4 / T(6) - T(4) * Q6 + Q3 / Q7)) / T(4);
  roots[3] = (-b + Q7 + complex_sqrt(T(4) * Q4 / T(6) - T(4) * Q6 + Q3 / Q7)) / T(4);
}

}

#endif
