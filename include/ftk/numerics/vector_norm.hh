#ifndef _FTK_NORM_H
#define _FTK_NORM_H

namespace ftk {

template <int n, typename T>
inline T vector_norm2(const T v[])
{
  T norm(0);
  for (int i=0; i<n; i++) 
    norm += v[i] * v[i];
  return norm;
}

template <typename T>
inline T vector_norm2_3(const T v[])
{
  return vector_norm2<3,T>(v);
}

template <typename T>
inline T vector_norm2_4(const T v[])
{
  return vector_norm2<3,T>(v);
}

}

#endif
