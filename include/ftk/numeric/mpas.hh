#ifndef _FTK_MPAS_HH
#define _FTK_MPAS_HH

#include <ftk/config.hh>
#include <ftk/numeric/cross_product.hh>

namespace ftk {

template <typename F>
bool point_in_mpas_cell(const int nverts, const F Xv[][3], const F x[3])
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

}

#endif
