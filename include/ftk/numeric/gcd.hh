#ifndef _FTK_GCD_HH
#define _FTK_GCD_HH

#include <ftk/ftk_config.hh>

namespace ftk {

template <typename I>
__device__ __host__ I gcd(I a, I b) {
  while (true) {
    if (a == I(0)) return b;
    b %= a;
    if (b == I(0)) return a;
    a %= b;
  }
}

} // namespace ftk

#endif
