#ifndef _FTK_VECTOR_ASSIGN_HH
#define _FTK_VECTOR_ASSIGN_HH

namespace ftk {

template <typename T>
inline void vector_assign3(T v[3], T x)
{
  v[0] = v[1] = v[2] = x;
}

template <typename T>
inline void vector_assign3(T v[3], T x, T y, T z)
{
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

template <typename T>
inline void vector_assign3(T v[3], const T w[3])
{
  v[0] = w[0];
  v[1] = w[1]; 
  v[2] = w[2];
}

// x = a*w
template <typename T>
inline void vector_assign_scalar_multiplication3(T x[3], T a, const T w[3])
{
  x[0] = a * w[0];
  x[1] = a * w[1];
  x[2] = a * w[2];
}

// x = w + v;
template <typename T>
inline void vector_assign_addition3(T x[3], const T w[3], const T v[3])
{
  x[0] = w[0] + v[0];
  x[1] = w[1] + v[1];
  x[2] = w[2] + v[2];
}

// x = w - v;
template <typename T>
inline void vector_assign_subtraction3(T x[3], const T w[3], const T v[3])
{
  x[0] = w[0] - v[0];
  x[1] = w[1] - v[1];
  x[2] = w[2] - v[2];
}

// x = a*w + b*v
template <typename T>
inline void vector_assign_weighted_addition3(T x[3], T a, const T w[3], T b, const T v[3])
{
  x[0] = a * w[0] + b * v[0];
  x[1] = a * w[1] + b * v[1];
  x[2] = a * w[2] + b * v[2];
}

}

#endif
