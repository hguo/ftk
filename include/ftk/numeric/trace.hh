#ifndef _FTK_TRACE_H
#define _FTK_TRACE_H

#include <ftk/ftk_config.hh>

namespace ftk {

template <typename T, int n>
__device__ __host__
inline T trace(T A[n][n])
{
  T tr(0);
  for (int i = 0; i < n; i ++)
    tr += A[i][i];
  return tr;
}

template <typename T>
__device__ __host__
inline T trace2(T A[2][2])
{
  return A[0][0] + A[1][1];
}

template <typename T>
__device__ __host__
inline T trace3(T m[3][3])
{
  return m[0][0] + m[1][1] + m[2][2];
}

}

#endif
