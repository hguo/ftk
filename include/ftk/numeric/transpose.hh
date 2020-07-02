#ifndef _FTK_TRANSPOSE_H
#define _FTK_TRANSPOSE_H

#include <ftk/ftk_config.hh>
#include <algorithm>
#include <utility>

namespace ftk {

template <typename T, int m, int n> // transpose a mxn matrix into a nxm matrix
__device__ __host__
void transpose(const T A[m][n], T B[n][m])
{
  for (int i = 0; i < m; i ++)
    for (int j = 0; j < n; j ++)
      B[j][i] = A[i][j];
}

template <typename T, int n> // transpose an nxn matrix in place
void transpose(T M[n][n])
{
  for (int i = 0; i < n; i ++)
    for (int j = i+1; j < n; j ++)
      std::swap(M[i][j], M[j][i]);
}

template <typename T>
__device__ __host__
void transpose2x2(const T A[2][2], T B[2][2])
{
  B[0][0] = A[0][0];
  B[0][1] = A[1][0];
  B[1][0] = A[0][1];
  B[1][1] = A[1][1];
}

template <typename T>
__device__ __host__
void transpose2x2(T m[4])
{
  std::swap(m[1], m[2]);
}

template <typename T>
__device__ __host__
void transpose3x3(T m[9]) 
{
  std::swap(m[1], m[3]);
  std::swap(m[2], m[6]);
  std::swap(m[5], m[7]);
}

template <typename T>
__device__ __host__
void transpose3x3(const T a[9], T b[9]) 
{
  b[0] = a[0];
  b[1] = a[3];
  b[2] = a[6];
  b[3] = a[1];
  b[4] = a[4];
  b[5] = a[7];
  b[6] = a[2];
  b[7] = a[5];
  b[8] = a[8];
}

template <typename T>
__device__ __host__
void transpose3x3(T A[3][3])
{
  std::swap(A[0][1], A[1][0]);
  std::swap(A[0][2], A[2][0]);
  std::swap(A[1][2], A[2][1]);
}

template <typename T>
__device__ __host__
void transpose3x3(const T A[3][3], T B[3][3])
{
  B[0][0] = A[0][0];
  B[0][1] = A[1][0];
  B[0][2] = A[2][0];
  B[1][0] = A[0][1];
  B[1][1] = A[1][1];
  B[1][2] = A[2][1];
  B[2][0] = A[0][2];
  B[2][1] = A[1][2];
  B[2][2] = A[2][2];
}

template <typename T>
void transpose4x4(T m[16]) 
{
  std::swap(m[1], m[4]);
  std::swap(m[2], m[8]);
  std::swap(m[3], m[12]);
  std::swap(m[6], m[9]);
  std::swap(m[7], m[13]);
  std::swap(m[11], m[14]);
}

template <typename T>
__device__ __host__
void transpose3x2(const T A[3][2], T B[2][3])
{
  B[0][0] = A[0][0]; 
  B[0][1] = A[1][0]; 
  B[0][2] = A[2][0];
  B[1][0] = A[0][1];
  B[1][1] = A[1][1];
  B[1][2] = A[2][1];
}

template <typename T>
__device__ __host__
void transpose2x3(const T A[2][3], T B[3][2])
{
  B[0][0] = A[0][0]; 
  B[0][1] = A[1][0];
  B[1][0] = A[0][1];
  B[1][1] = A[1][1];
  B[2][0] = A[0][2];
  B[2][1] = A[1][2];
}

template <typename T>
__device__ __host__
void transpose4x3(const T A[4][3], T B[3][4])
{
  B[0][0] = A[0][0]; 
  B[0][1] = A[1][0]; 
  B[0][2] = A[2][0];
  B[0][3] = A[3][0];

  B[1][0] = A[0][1];
  B[1][1] = A[1][1];
  B[1][2] = A[2][1];
  B[1][3] = A[3][1];

  B[2][0] = A[0][2];
  B[2][1] = A[1][2];
  B[2][2] = A[2][2];
  B[2][3] = A[3][2];
}

}

#endif
