#ifndef _FTK_NORMALIZE_H
#define _FTK_NORMALIZE_H

#include <ftk/numeric/vector_norm.hh>

namespace ftk {

template <int n, typename T>
inline T vector_normalization2(T v[])
{
  T norm = vector_2norm<n, T>(v);
  for (int i=0; i<n; i++) 
    v[i] /= norm;
  return norm;
}

template <typename T>
inline T vector_normalization2_2(T v[2])
{
  return vector_normalization2<2,T>(v);
}

template <typename T>
inline void vector_normalization2_2(const T v[2], T w[2])
{
  w[0] = v[0]; w[1] = v[1];
  vector_normalization2_2(w);
}

template <typename T>
inline T vector_normalization2_3(T v[])
{
  return vector_normalization2<3,T>(v);
}

template <typename T>
inline T vector_normalization2_3(const T v[3], T w[3])
{
  w[0] = v[0]; w[1] = v[1]; w[2] = v[2];
  return vector_normalization2_3(w);
}

template <typename T>
inline T vector_normalization2_4(T v[])
{
  return vector_normalization2<4,T>(v);
}

}

#endif
