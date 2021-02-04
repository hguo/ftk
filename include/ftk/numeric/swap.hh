#ifndef _FTK_SWAP_HH
#define _FTK_SWAP_HH

#include <ftk/config.hh>

template <typename T>
__device__ __host__
inline void swap_helper(T& a, T& b)
{
  T c(a); a = b; b = c;
}

#endif
