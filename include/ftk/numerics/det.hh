#ifndef _FTK_INVMAT_H
#define _FTK_INVMAT_H

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

}

#endif
