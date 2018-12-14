#ifndef _FTK_POLYNOMIAL_HH
#define _FTK_POLYNOMIAL_HH

namespace ftk {

template <typename T>
void polynomial_copy(const T P[], int m, T Q[])
{
  for (int i = 0; i < m+1; i ++) 
    Q[i] = P[i];
}

template <typename T>
void polynomial_add(const T P[], int m, const T Q[], int n, T R[]) 
{
  if (m >= n)
    for (int i = 0; i < m+1; i ++)
      if (i < n) R[i] = P[i] + Q[i];
      else R[i] = P[i];
  else 
    polynomial_add(Q, n, P, m, R);
}

template <typename T>
void polynomial_add_in_place(T P[], int m, const T Q[], int n) // m >= n
{
  for (int i = 0; i < n+1; i ++)
    P[i] += Q[i];
}

template <typename T>
void polynomial_multiplication(const T P[], int m, const T Q[], int n, T R[])
{
  for (int i = 0; i < m+n+1; i ++) 
    R[i] = T(0);

  for (int i = 0; i < m+1; i ++) 
    for (int j = 0; j < n+1; j ++) 
      R[i+j] += P[i] * Q[j];
}

template <typename T>
void polynomial_scalar_multiplication(T P[], int m, T scalar)
{
  for (int i = 0; i < m+1; i ++)
    P[i] *= scalar;
}

template <typename T>
void polynomial_scalar_multiplication(const T P[], int m, T scalar, const T R[])
{
  for (int i = 0; i < m+1; i ++) 
    R[i] = P[i] * scalar;
}

}

#endif
