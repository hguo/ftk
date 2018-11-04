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

}

#endif
