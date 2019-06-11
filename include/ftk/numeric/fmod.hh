#ifndef _FTK_FMOD_H
#define _FTK_FMOD_H

#include <cmath>

namespace ftk {

template <typename T>
inline static T fmod1(T x, T y)
{
  T z = fmod(x, y);
  if (z<0) z += y;
  return z;
}

template <typename T>
inline static T mod2pi(T x)
{
  T y = fmod(x, 2*M_PI); 
  if (y<0) y+= 2*M_PI;
  return y; 
}

template <typename T>
inline static T mod2pi1(T x)
{
  return mod2pi(x + M_PI) - M_PI;
}

}

#endif
