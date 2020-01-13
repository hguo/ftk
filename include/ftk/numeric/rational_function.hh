#ifndef _FTK_RATIONAL_FUNCTION_HH
#define _FTK_RATIONAL_FUNCTION_HH

#include <ftk/numeric/polynomial.hh>

namespace ftk {

template <typename T>
T evaluate_rational(const T P[], const T Q[], int n, T x, const T epsilon = std::numeric_limits<T>::epsilon())
{
  // fprintf(stderr, "P=%f, Q=%f\n", 
  //     polynomial_evaluate(P, n, x), 
  //     polynomial_evaluate(Q, n, x));

  if (std::abs(polynomial_evaluate(Q, n, x)) < epsilon) {
    if (n > 0 && std::abs(polynomial_evaluate(P, n, x)) < epsilon) { // P=0, Q=0, use L'Hospital's rule
      T dPdx[n], dQdx[n];
      polynomial_derivative(P, n, dPdx);
      polynomial_derivative(Q, n, dQdx);
      return evaluate_rational(dPdx, dQdx, n-1, x, epsilon); // dPdx/dQdx
    } else 
      return std::nan(""); // P!=0, Q=0
  } else return polynomial_evaluate(P, n, x) / polynomial_evaluate(Q, n, x);
}

template <typename T>
T evaluate_rational_infinity(const T P[], const T Q[], int n, const T epsilon = std::numeric_limits<T>::epsilon())
{
  if (std::abs(Q[n]) < epsilon) {
    if (n > 0 && std::abs(P[n]) < epsilon) 
      return evaluate_rational_infinity(P, Q, n-1, epsilon);
    else return std::nan("");
  } else return P[n] / Q[n];
}

#if 0
template <typename T>
T evaluate_cubic_rational(const T P[4], const T Q[4], const T x, const T epsilon = 1e-9)
{
  if (std::abs(polynomial_evaluate(Q, 3, x)) < epsilon) {
    if (std::abs(polynomial_evaluate(P, 3, x)) < epsilon) {

    }
  }

#if 0
  if (polynomial_equals_to_zero(Q, 3, epsilon)) return std::nan(""); // the denominator constantly equals to zero
  else if (polynomial_equals_to_zero(P, 3, epsilon)) return T(0); // the numerator constantly equals to zero but the denominator is not.
  else if (std::abs(polynomial_evaluate(Q, 3, x)) < epsilon) 
#endif
}

template <typename T>
T evaluate_linear_rational_infinity(const T P[2], const T Q[2], const T epsilon = 1e-9)
{
  if (std::abs(Q[1]) < epsilon) {
    if (std::abs(P[1]) < epsilon) 
      return P[0] / Q[0];
    else 
      return std::nan("");
  } else return 
    P[1] / Q[1];
}

template <typename T>
T evaluate_quadratic_rational_infinity(const T P[3], const T Q[3], const T epsilon = 1e-9)
{
  if (std::abs(Q[2]) < epsilon) {
    if (std::abs(P[2]) < epsilon)
      return evaluate_linear_rational_infinity(P, Q, epsilon);
    else 
      return std::nan("");
  } else 
    return P[2] / Q[2];
}

template <typename T>
T evaluate_cubic_rational_infinity(const T P[4], const T Q[4], const T epsilon = 1e-9)
{
  if (std::abs(Q[3]) < epsilon) {
    if (std::abs(P[3]) < epsilon)
      return evaluate_quadratic_rational_infinity(P, Q, epsilon);
    else 
      return std::nan("");
  } else 
    return P[3] / Q[3];
}

template <typename T>
T evaluate_cubic_rational_zero(const T P[4], const T Q[4])
{
  if (std::abs(Q[0]) < epsilon) {
    if (std::abs(P[0]) < epsilon)
      return evaluate_quadratic_rational_infinity(P+1, Q+1, epsilon);
    else return std::nan("");
  } else
    return P[0] / Q[0];
}
#endif

}

#endif
