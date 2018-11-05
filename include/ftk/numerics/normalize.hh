#ifndef _FTK_NORMALIZE_H
#define _FTK_NORMALIZE_H

#include <ftk/numerics/norm.hh>

namespace ftk {

template <int n, typename T>
inline void normalize2_vector(T v[])
{
  T norm = norm2_vector<n, T>(v);
  for (int i=0; i<n; i++) 
    v[i] /= norm;
}

template <typename T>
inline void normalize2_3(T v[])
{
  normalize2_vector<3,T>(v);
}

template <typename T>
inline void normalize2_4(T v[])
{
  normalize2_vector<3,T>(v);
}

}

#endif
