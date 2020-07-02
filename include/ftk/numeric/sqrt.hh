#ifndef _FTK_SQRT_HH
#define _FTK_SQRT_HH

#include <cmath>
#include <complex>

namespace ftk {

template <typename T>
__device__ __host__
static std::complex<T> complex_sqrt(const std::complex<T> z)
{
  return pow(z, T(1)/T(2));
}

}

#endif
