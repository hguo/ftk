#ifndef _FTK_LINEAR_SOLVE1_H
#define _FTK_LINEAR_SOLVE1_H

namespace ftk {

template <typename T>
inline int solve_linear_real1(const T P[2], T x[1], const T epsilon = std::numeric_limits<T>::epsilon())
{
  if (std::abs(P[1]) < epsilon || std::isinf(P[1]) || std::isnan(P[1])) return 0;
  else {
    if (std::isinf(P[0]) || std::isnan(P[0])) return 0; 
    else {
      x[0] = - P[0] / P[1];
      return 1;
    }
  }
}

template <typename T>
inline T solve_least_square2x1(const T a[2], const T b[2], T &x) // ax=b
{
  const T inner_product_a = inner_product2(a, a),
          inner_product_ab = inner_product2(a, b);
  const T inv_inner_product_a = T(1) / inner_product_a;
  x = inner_product_ab * inv_inner_product_a;
  return inv_inner_product_a;
}

template <typename T>
inline T solve_least_square3x1(const T a[3], const T b[3], T &x) // ax=b
{
  const T inner_product_a = inner_product3(a, a), 
          inner_product_ab = inner_product3(a, b);
  const T inv_inner_product_a = T(1) / inner_product_a;
  x = inner_product_ab * inv_inner_product_a;
  return inv_inner_product_a; // cond
}

}

#endif
