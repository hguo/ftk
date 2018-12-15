#ifndef _FTK_POLYNOMIAL_HH
#define _FTK_POLYNOMIAL_HH

#include <cmath>

namespace ftk {

template <typename T>
T polynomial_evaluate(const T P[], int m, T x)
{
  T y(0);
  for (int i = 0; i <= m; i ++) 
    y += P[i] * std::pow(x, i);
  return y;
}

template <typename T>
void polynomial_derivative(const T P[], int m, T R[])
{
  for (int i = 0; i < m; i ++) 
    R[i] = P[i+1] * (i+1);
}

template <typename T>
T polynomial_derivative_evaluate(const T P[], int m, T x)
{
  T y(0);
  for (int i = 0; i < m; i++)
    y += P[i+1] * (i+1) * std::pow(x, i);
  return y;
}

template <typename T>
void polynomial_copy(const T P[], int m, T Q[])
{
  for (int i = 0; i <= m; i ++) 
    Q[i] = P[i];
}

template <typename T>
void polynomial_add(const T P[], int m, const T Q[], int n, T R[]) 
{
  if (m >= n)
    for (int i = 0; i <= m; i ++)
      if (i <= n) R[i] = P[i] + Q[i];
      else R[i] = P[i];
  else 
    polynomial_add(Q, n, P, m, R);
}

template <typename T>
void polynomial_add_in_place(T P[], int m, const T Q[], int n) // m >= n
{
  for (int i = 0; i <= n; i ++)
    P[i] += Q[i];
}

template <typename T>
void polynomial_multiplication(const T P[], int m, const T Q[], int n, T R[])
{
  for (int i = 0; i <= m+n; i ++) 
    R[i] = T(0);

  for (int i = 0; i <= m; i ++) 
    for (int j = 0; j <= n; j ++) 
      R[i+j] += P[i] * Q[j];
}

template <typename T>
void polynomial_scalar_multiplication(T P[], int m, T scalar)
{
  for (int i = 0; i <= m; i ++)
    P[i] *= scalar;
}

template <typename T>
void polynomial_scalar_multiplication(const T P[], int m, T scalar, const T R[])
{
  for (int i = 0; i <= m; i ++) 
    R[i] = P[i] * scalar;
}

}

#endif
