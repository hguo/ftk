#ifndef _FTK_NORM_H
#define _FTK_NORM_H

namespace ftk {

template <int n, typename ValueType>
inline ValueType vecnorm2(const ValueType v[])
{
  ValueType norm(0);
  for (int i=0; i<n; i++) 
    norm = v[i] * v[i];
  return norm;
}

}

#endif
