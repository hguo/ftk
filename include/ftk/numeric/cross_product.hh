#ifndef _FTK_CROSS_PRODUCT_H
#define _FTK_CROSS_PRODUCT_H

namespace ftk {

template <typename T>
static inline T cross_product2(const T A[2], const T B[2])
{
  return A[0]*B[1] - A[1]*B[0];
}

template <typename T>
static inline void cross_product(const T A[3], const T B[3], T C[3])
{
  C[0] = A[1]*B[2] - A[2]*B[1]; 
  C[1] = A[2]*B[0] - A[0]*B[2]; 
  C[2] = A[0]*B[1] - A[1]*B[0];
}

}

#endif
