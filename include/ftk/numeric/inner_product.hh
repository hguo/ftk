#ifndef _FTK_INNER_PRODUCT_H
#define _FTK_INNER_PRODUCT_H

namespace ftk {

template <typename T>
static inline T inner_product2(const T A[2], const T B[2])
{
  return A[0]*B[0] + A[1]*B[1];
}

template <typename T>
static inline T inner_product3(const T A[3], const T B[3])
{
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

}

#endif
