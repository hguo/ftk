#ifndef _FTK_CLAMP_HH
#define _FTK_CLAMP_HH

#include <ftk/config.hh>

namespace ftk {

#ifndef __CUDACC__
  using std::min;
  using std::max;
  using std::isnan;
  using std::isinf;
#endif

template <typename T>
__device__ __host__
inline T clamp(T x, T a, T b)
{
  return min(max(a, x), b);
}

template <int n, typename T>
__device__ __host__
inline void clamp_barycentric(T x[])
{
  T sum = 0.0;
  for (int i = 0; i < n; i ++) {
    x[i] = clamp<T>(x[i], 0.0, 1.0);
    sum += x[i];
  }
  for (int i = 0; i < n; i ++)
    x[i] /= sum;

  if (isnan(x[0]) || isinf(x[0]))
    for (int i = 0; i < n; i ++)
      x[i] = 1.0 / n;
}

}

#endif
