#ifndef _FTK_NORM_H
#define _FTK_NORM_H

namespace ftk {

template <int n, typename T>
inline T norm2_vector(const T v[])
{
  T norm(0);
  for (int i=0; i<n; i++) 
    norm = v[i] * v[i];
  return norm;
}

template <typename T>
inline T norm2_3(const T v[])
{
  return norm2_vector<3,T>(v);
}

template <typename T>
inline T norm2_4(const T v[])
{
  return norm2_vector<3,T>(v);
}

}

#endif
