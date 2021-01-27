#ifndef _FTK_DET_H
#define _FTK_DET_H

#include <ftk/config.hh>
#include <iostream>

namespace ftk {

template <typename T>
__device__ __host__
inline T det2(const T A[2][2])
{
  return A[0][0] * A[1][1] - A[1][0] * A[0][1];
}

template <class T>
__device__ __host__
inline T det3(const T m[3][3])
{
  return 
      m[0][0] * (m[1][1]*m[2][2] - m[1][2]*m[2][1])
    - m[0][1] * (m[1][0]*m[2][2] - m[1][2]*m[2][0])
    + m[0][2] * (m[1][0]*m[2][1] - m[1][1]*m[2][0]);
}

template <class T>
__device__ __host__
inline T det4(const T m[4][4])
{
  const T 
    d2233 = m[2][2] * m[3][3] - m[2][3] * m[3][2],
    d2133 = m[2][1] * m[3][3] - m[2][3] * m[3][1],
    d2132 = m[2][1] * m[3][2] - m[2][2] * m[3][1],
    d2033 = m[2][0] * m[3][3] - m[2][3] * m[3][0],
    d2032 = m[2][0] * m[3][2] - m[2][2] * m[3][0],
    d2031 = m[2][0] * m[3][1] - m[2][1] * m[3][0];

  return
      m[0][0] * (
            m[1][1] * d2233 
          - m[1][2] * d2133 
          + m[1][3] * d2132)
    - m[0][1] * (
            m[1][0] * d2233
          - m[1][2] * d2033
          + m[1][3] * d2032)
    + m[0][2] * (
            m[1][0] * d2133
          - m[1][1] * d2033
          + m[1][3] * d2031)
    - m[0][3] * (
            m[1][0] * d2132
          - m[1][1] * d2032
          + m[1][2] * d2031);
}

template <class T>
__device__ __host__
inline T det2(T m00, T m01, T m10, T m11)
{
  return m00*m11 - m10*m01;
}

}

#endif
