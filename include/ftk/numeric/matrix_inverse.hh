#ifndef _FTK_MATRIX_INVERSE_H
#define _FTK_MATRIX_INVERSE_H

#include <ftk/ftk_config.hh>

namespace ftk {

template <class T>
__device__ __host__
inline T matrix_inverse2x2(const T m[2][2], T inv[2][2])
{
  const T det = m[0][0]*m[1][1] - m[0][1]*m[1][0];
  const T invdet = T(1) / det;
  
  inv[0][0] = m[1][1] * invdet;
  inv[0][1] = -m[0][1] * invdet;
  inv[1][0] = -m[1][0] * invdet;
  inv[1][1] = m[0][0] * invdet;

  return det;
}

template <class T>
__device__ __host__
inline T matrix_inverse3x3(const T m[3][3], T inv[3][3]) // returns the determinant
{
  inv[0][0] =   m[1][1]*m[2][2] - m[1][2]*m[2][1];
  inv[0][1] = - m[0][1]*m[2][2] + m[0][2]*m[2][1];
  inv[0][2] =   m[0][1]*m[1][2] - m[0][2]*m[1][1];
  inv[1][0] = - m[1][0]*m[2][2] + m[1][2]*m[2][0];
  inv[1][1] =   m[0][0]*m[2][2] - m[0][2]*m[2][0];
  inv[1][2] = - m[0][0]*m[1][2] + m[0][2]*m[1][0];
  inv[2][0] =   m[1][0]*m[2][1] - m[1][1]*m[2][0];
  inv[2][1] = - m[0][0]*m[2][1] + m[0][1]*m[2][0];
  inv[2][2] =   m[0][0]*m[1][1] - m[0][1]*m[1][0];
  
  T det = m[0][0]*inv[0][0] + m[0][1]*inv[1][0] + m[0][2]*inv[2][0];
  T invdet = T(1) / det;

  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      inv[i][j] = inv[i][j] * invdet;

  return det;
}

template <class T>
__device__ __host__
inline T matrix_inverse4(const T m[16], T inv[16]) // returns det
{
  inv[0] =   m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15]
           + m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
  inv[4] =  -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15]
           - m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
  inv[8] =   m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15]
           + m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
  inv[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14]
           - m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
  inv[1] =  -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15]
           - m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
  inv[5] =   m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15]
           + m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
  inv[9] =  -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15]
           - m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
  inv[13] =  m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14]
           + m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
  inv[2] =   m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15]
           + m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
  inv[6] =  -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15]
           - m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
  inv[10] =  m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15]
           + m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
  inv[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14]
           - m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];
  inv[3] =  -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11]
           - m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
  inv[7] =   m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11]
           + m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
  inv[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11]
           - m[4]*m[3]*m[9] - m[8]*m[1]*m[7] + m[8]*m[3]*m[5];
  inv[15] =  m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10]
           + m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];

  T det = m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12];
  T invdet = T(1) / det;

  for (int i = 0; i < 16; i++)
    inv[i] = inv[i] * invdet;

  return det;
}

}

#endif
