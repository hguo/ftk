#ifndef _FTK_SYMMETRIC_MATRIX_HH
#define _FTK_SYMMETRIC_MATRIX_HH

#include <ftk/ftk_config.hh>

namespace ftk {

template <typename T>
__device__ __host__
inline void make_symmetric2x2(T m[2][2])
{
  T a = T(0.5) * (m[0][1] + m[1][0]);
  m[0][1] = m[1][0] = a;
}

}

#endif
