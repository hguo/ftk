#ifndef _FTK_DET_H
#define _FTK_DET_H

#include <iostream>

namespace ftk {

template <class ValueType>
inline ValueType det2(const ValueType m[4])
{
  return m[0]*m[3] - m[1]*m[2];
}

template <class ValueType>
inline ValueType det3(const ValueType m[9])
{
  return m[0] * (m[4]*m[8] - m[5]*m[7])
    + m[1] * (-m[3]*m[8] + m[5]*m[6])
    + m[2] * (m[3]*m[7] - m[4]*m[6]);
}

template <class ValueType>
inline ValueType det3(const ValueType m[3][3]) // untested
{
  return m[0][0] * (m[1][1]*m[2][2] - m[1][2]*m[2][1])
    + m[0][1] * (-m[1][0]*m[2][2] + m[1][2]*m[2][0])
    + m[0][2] * (m[1][0]*m[2][1] - m[1][1]*m[2][0]);
}

template <class ValueType>
inline ValueType det4(const ValueType m[16])
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
inline T det2(T m0, T m1, T m2, T m3)
{
  T m[] = {m0, m1, m2, m3};
  return det2(m);
}

template <class T>
inline T det3(T m0, T m1, T m2, T m3, T m4, T m5, T m6, T m7, T m8)
{
  T m[] = {m0, m1, m2, m3, m4, m5, m6, m7, m8};
  return det3(m);
}

}

#endif
