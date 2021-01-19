#ifndef _FTK_XGC_INTERPOLANT_HH
#define _FTK_XGC_INTERPOLANT_HH

namespace ftk {

template <typename I=int, typename F=double>
struct xgc_interpolant_t { // poloidal interpolants
  I tri0[3], tri1[3];
  F mu0[3], mu1[3];
};

} // namespace ftk

#endif
