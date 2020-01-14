#ifndef _FTK_ORIENTATION_HH
#define _FTK_ORIENTATION_HH

#include <ftk/ftk_config.hh>
#include <ftk/sign_det.hh>

namespace ftk {

// unrobust orientation test
template <typename T=long long>
inline int orientation2(const T x0[2], const T x1[2], const T x2[2])
{
  T X[3][3] = {
    {x0[0], x0[1], T(1)}, 
    {x1[0], x1[1], T(1)},
    {x2[0], x2[1], T(1)}
  };
  return sign_det3(X);
}

} // namespace ftk

#endif
