#ifndef _ACCESS_H
#define _ACCESS_H

namespace ftk {

template <class IdType, class NumberType, class ArrayType>
NumberType& access2D(ArrayType &arr, IdType W, IdType H, IdType i, IdType j)
{
  return arr[i + W*j];
}

}

#if 0
inline float texel2D(const float *ptr, const int *sz, int x, int y);
inline float texel3D(const float *ptr, const int *sz, int x, int y, int z);
inline float texel4D(const float *ptr, const int *sz, int x, int y, int z, int t);

//////////
inline float texel2D(const float* p, const int* sz, int x, int y)
{
  return p[x + sz[0]*y];
}

inline float texel3D(const float* p, const int* sz, int x, int y, int z)
{
  return p[x + sz[0]*(y + sz[1]*z)];
}

inline float texel4D(const float* p, const int* sz, int x, int y, int z, int t)
{
  return p[x + sz[0]*(y + sz[1]*(z + sz[2]*t))];
}
#endif

#endif
