#ifndef _FTK_ISNAN_H
#define _FTK_ISNAN_H

#include <cmath>

namespace ftk {

template <typename T>
inline bool isnan3x3(const T m[])
{
  for (int i=0; i<9; i++) 
    if (std::isnan(m[i])) return true;
  return false;
}

template <typename T>
inline bool isnan3x3(const T m[3][3])
{
  for (int i=0; i<3; i++) 
    for (int j=0; j<3; j++) 
      if (std::isnan(m[i][j])) return true;
  return false;
}

}

#endif
