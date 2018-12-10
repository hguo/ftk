#ifndef _FTK_NORMALIZE_H
#define _FTK_NORMALIZE_H

#include <ftk/numeric/vector_norm.hh>

namespace ftk {

template <int n, typename T>
inline void vector_normalization2(T v[])
{
  T norm = vector_2norm<n, T>(v);
  for (int i=0; i<n; i++) 
    v[i] /= norm;
}

template <typename T>
inline void vector_normalization2_3(T v[])
{
  vector_normalization2<3,T>(v);
}

template <typename T>
inline void vector_normalization2_4(T v[])
{
  vector_normalization2<4,T>(v);
}

}

#endif
