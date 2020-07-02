#ifndef _FTK_SIGN_HH
#define _FTK_SIGN_HH

#include <ftk/ftk_config.hh>

namespace ftk {

template <typename T>
__device__ __host__
inline int sign(T x)
{
  // https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
  return (T(0) < x) - (x < T(0));
}

} // namespace ftk

#endif
