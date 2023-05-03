#ifndef _FTK_WACHPRESS_INTERPOLATION_HH
#define _FTK_WACHPRESS_INTERPOLATION_HH

#include <ftk/config.hh>

namespace ftk {

template <typename T=double, int nd=3, int nvars=3>
__device__ __host__ 
inline bool wachpress_coords(
    const int nverts,
    const T X[][nd], const T V[][nvars], // inputs
    const T x[nd], // cartisian coordinates
    T mu[nvars]) // wachpress coordinates
{
  return false;
}

template <typename T=double, int nd=3, int nvars=3>
__device__ __host__ 
inline bool wachpress_interpolation(
    const int nverts,
    const T X[][nd], const T V[][nvars], // inputs
    const T x[nd], // coordinates
    T f[nvars]) // outputs
{
  return false;
}

}

#endif
