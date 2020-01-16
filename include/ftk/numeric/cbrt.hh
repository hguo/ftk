#ifndef _FTK_CBRT_HH
#define _FTK_CBRT_HH

#include <cmath>
#include <complex>

namespace ftk {

template <typename T>
__device__ __host__
static std::complex<T> complex_cbrt(const std::complex<T> z)
{
  return pow(z, T(1)/T(3));
}

}

#endif
