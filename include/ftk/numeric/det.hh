#ifndef _FTK_DET_H
#define _FTK_DET_H

#include <ftk/ftk_config.hh>
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
inline T det3(const T m[3][3]) // untested
{
  return m[0][0] * (m[1][1]*m[2][2] - m[1][2]*m[2][1])
    + m[0][1] * (-m[1][0]*m[2][2] + m[1][2]*m[2][0])
    + m[0][2] * (m[1][0]*m[2][1] - m[1][1]*m[2][0]);
}

template <class T>
__device__ __host__
inline T det4(const T m[16]) // legacy code
{
  return 
      m[1] * m[11] * m[14] * m[4] 
    - m[1] * m[10] * m[15] * m[4] 
    - m[11] * m[13] * m[2] * m[4] 
    + m[10] * m[13] * m[3] * m[4] 
    - m[0] * m[11] * m[14] * m[5]
    + m[0] * m[10] * m[15] * m[5] 
    + m[11] * m[12] * m[2] * m[5] 
    - m[10] * m[12] * m[3] * m[5] 
    - m[1] * m[11] * m[12] * m[6] 
    + m[0] * m[11] * m[13] * m[6]
    + m[1] * m[10] * m[12] * m[7] 
    - m[0] * m[10] * m[13] * m[7] 
    - m[15] * m[2] * m[5] * m[8] 
    + m[14] * m[3] * m[5] * m[8]
    + m[1] * m[15] * m[6] * m[8]
    - m[13] * m[3] * m[6] * m[8]
    - m[1] * m[14] * m[7] * m[8] 
    + m[13] * m[2] * m[7] * m[8]
    + m[15] * m[2] * m[4] * m[9]
    - m[14] * m[3] * m[4] * m[9] 
    - m[0] * m[15] * m[6] * m[9]
    + m[12] * m[3] * m[6] * m[9] 
    + m[0] * m[14] * m[7] * m[9]
    - m[12] * m[2] * m[7] * m[9];
}

template <class T>
__device__ __host__
inline T det2(T m00, T m01, T m10, T m11)
{
  return m00*m11 - m10*m01;
}

}

#endif
