#ifndef _FTK_MPAS_HH
#define _FTK_MPAS_HH

#include <ftk/config.hh>
#include <ftk/numeric/cross_product.hh>

namespace ftk {

template <typename F=double>
__device__ __host__ 
inline bool point_in_mpas_cell(const int nverts, const F Xv[][3], const F x[3])
{
  F n[3];
  for (int i = 0; i < nverts; i ++) {
    const F *x0 = Xv[i], *x1 = Xv[(i+1)%nverts];
    cross_product(x0, x1, n);

    if (vector_dot_product3(n, x) < 0) // on the negative side
      return false;
  }

  return true;
}

template <typename F=double>
__device__ __host__
inline int locate_mpas_layer_brute_force(const int nlayers, const F zTop[], const F z)
{
  int i0 = 0, i1 = nlayers;
  // I ib = locate_layer_bisection(nlayers, zTop, z);

  // fprintf(stderr, "z=%f, ib=%d\n", z, ib);
  // for (int i = 0; i < nlayers; i ++)
  //   fprintf(stderr, "i=%d, ztop=%f\n", i, zTop[i]);

  for (int i = 0; i < nlayers; i ++)
    if (z <= zTop[i] && z > zTop[i+1])
      return i;

  return -1;
}

template <typename F=double>
__device__ __host__
inline int locate_mpas_layer_bisection(const int nlayers, const F zTop[], const F z)
{
  int i0 = 0, i1 = nlayers-1;
  // fprintf(stderr, "z=%f, zTop_i0=%f, zTop_i1=%f\n", z, zTop[i0], zTop[i1]);

  if (z > zTop[i0] || z <= zTop[i1])
    return -1;

#if 0
  fprintf(stderr, "--z=%f\n", z);
  for (int i = 0; i < nlayers; i ++)
    fprintf(stderr, "--i=%d, ztop=%f\n", i, zTop[i]);
  // exit(1);
#endif

  while (1) {
    if (z <= zTop[i0] && z > zTop[i0+1])
      return i0;
    
    int ic = (i0 + i1) / 2;
    // fprintf(stderr, "ic=%d\n", ic);
    if (z - zTop[ic] <= 0)
      i0 = ic;
    else 
      i1 = ic;
  }

  return -1;
}

}

#endif
