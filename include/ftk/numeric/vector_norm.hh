#ifndef _FTK_NORM_H
#define _FTK_NORM_H

namespace ftk {

template <int n, typename T>
inline T vector_2norm(const T v[])
{
  T norm(0);
  for (int i=0; i<n; i++) 
    norm += v[i] * v[i];
  return sqrt(norm);
}

template <typename T>
inline T vector_2norm_3(const T v[])
{
  return vector_2norm<3,T>(v);
}

template <typename T>
inline T vector_2norm_4(const T v[])
{
  return vector_2norm<4,T>(v);
}

}

#endif
