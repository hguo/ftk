#ifndef _FTK_RATIONAL_HH
#define _FTK_RATIONAL_HH

#include <polynomial>

template <typename T>
T rational_evaluate(const T P[], int m, const T Q, int n, T x)
{
  return polynomial_evaluate(P, m, x) / polynomial_evaluate(Q, n, x);
}

#endif
