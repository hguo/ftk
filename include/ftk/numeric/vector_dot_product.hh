#ifndef _FTK_VECTOR_DOT_PRODUCT_HH
#define _FTK_VECTOR_DOT_PRODUCT_HH

namespace ftk {

template <typename T>
__device__ __host__
inline T vector_dot_product2(const T A[2], const T B[2])
{
  return A[0]*B[0] + A[1]*B[1];
}

template <typename T>
__device__ __host__
inline T vector_dot_product3(const T A[3], const T B[3])
{
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

}

#endif
