#ifndef _FTK_MATRIX_ADDITION_H
#define _FTK_MATRIX_ADDITION_H

#include <ftk/ftk_config.hh>

namespace ftk {

template <class T, int m, int n>
void matrix_addition(const T A[m][n], const T B[m][n], T C[m][n])
{
  for (int i = 0; i < m; i ++) 
    for (int j = 0; j < n; j ++)
      C[i][j] = A[i][j] + B[i][j];
}

template <class T, int m, int n>
void matrix_subtraction(const T A[m][n], const T B[m][n], T C[m][n])
{
  for (int i = 0; i < m; i ++) 
    for (int j = 0; j < n; j ++)
      C[i][j] = A[i][j] - B[i][j];
}

template <class T>
void matrix_addition3x3(const T A[3][3], const T B[3][3], T C[3][3])
{
  matrix_addition<T, 3, 3>(A, B, C);
}

template <class T>
void matrix_subtraction3x3(const T A[3][3], const T B[3][3], T C[3][3])
{
  matrix_subtraction<T, 3, 3>(A, B, C);
}

}

#endif
