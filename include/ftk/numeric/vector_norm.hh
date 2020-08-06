#ifndef _FTK_VECTOR_NORM_H
#define _FTK_VECTOR_NORM_H

#include <cmath>

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
inline T vector_2norm_2(const T v[])
{
  return std::sqrt(v[0]*v[0] + v[1]*v[1]);
}

template <typename T>
inline T vector_2norm_3(const T v[])
{
  return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  // return vector_2norm<3,T>(v);
}

template <typename T>
inline T vector_2norm_4(const T v[])
{
  return vector_2norm<4,T>(v);
}

template <typename T>
inline T vector_dist_2norm_2(const T v[], const T w[])
{
  T u[2] = {v[0] - w[0], v[1] - w[1]};
  return vector_2norm_2(u);
}

template <typename T>
inline T vector_dist_2norm_3(const T v[], const T w[])
{
  T u[3] = {v[0] - w[0], v[1] - w[1], v[2] - w[2]};
  return vector_2norm_3(u);
}

}

#endif
