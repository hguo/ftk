#ifndef _FTK_TRACE_H
#define _FTK_TRACE_H

namespace ftk {

template <typename T>
inline T trace2(T A[2][2])
{
  return A[0][0] + A[1][1];
}

template <typename T>
inline T trace3(T m[3][3])
{
  return m[0][0] + m[1][1] + m[2][2];
}

}

#endif
