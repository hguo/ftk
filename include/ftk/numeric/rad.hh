#ifndef _FTK_RAD_HH
#define _FTK_RAD_HH

#include <ftk/config.hh>

namespace ftk {

template <typename T=double>
__device__ __host__
inline T deg2rad(const T x)
{
  return x * M_PI / 180.0;
}

template <typename T=double>
__device__ __host__
inline T rad2deg(const T x)
{
  return x / M_PI * 180.0;
}

}

#endif
