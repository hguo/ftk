#ifndef _FTK_TRACE_H
#define _FTK_TRACE_H

namespace ftk {

template <typename ValueType>
inline ValueType trace3(ValueType m[]) 
{
  return m[0] + m[4] + m[8];
}

template <typename ValueType>
inline ValueType trace3(ValueType m[3][3]) // untested
{
  return m[0][0] + m[1][1] + m[2][2];
}

template <typename ValueType>
inline ValueType trace4(ValueType m[]) 
{
  return m[0] + m[5] + m[10] + m[15];
}

}

#endif
