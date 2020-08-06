#ifndef _FTK_LINEAR_SOLVE1_H
#define _FTK_LINEAR_SOLVE1_H

#include <cmath>
#include <limits>
#include <ftk/numeric/sign.hh>
#include <ftk/numeric/inner_product.hh>
#include <ftk/numeric/cross_product.hh>

namespace ftk {

enum {
  LINEAR_EQUATION_NO_SOLUTION = 0, 
  LINEAR_EQUATION_UNIQUE_SOLUTION = 1,
  LINEAR_EQUATION_INFINITE_SOLUTIONS = 
    std::numeric_limits<int>::max()
};

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
inline int number_roots_linear1(const T a, const T b) // ax=b, returns unique/infinite/no solution
{
  if (a == T(0)) {
    if (b == T(0)) return LINEAR_EQUATION_INFINITE_SOLUTIONS;
    else return LINEAR_EQUATION_NO_SOLUTION;
  } 
  else return LINEAR_EQUATION_UNIQUE_SOLUTION;
}

#if 0
template <typename T>
inline int number_roots_linear1(int n, const T a[], const T b[]) // check if n ax=b equations have a unique common root
{
  // if any of the equations has no solution, returns no solution.
  // if each of the equations has infinite solutions, return infinite solutions.
  // if each equation has the same solution, return unique solution; otherwise return no solution.
#if 0 // TODO
  bool all_infinite_solutions = true;
  for (int i = 0; i < n; i++) {
    const int ns = number_roots_linear1(a[i], b[i]);
    if (ns == LINEAR_EQUATION_NO_SOLUTION) return LINEAR_EQUATION_NO_SOLUTION;
    else if (ns == LINEAR_EQUATION_UNIQUE_SOLUTION) all_infinite_solutions = false;

  }
#endif
}
#endif

template <typename T>
inline T solve_least_square2x1(const T a[2], const T b[2], T &x) // ax=b
{
  const T inner_product_aa = inner_product2(a, a),
          inner_product_ab = inner_product2(a, b);
  const T inv_inner_product_aa = T(1) / inner_product_aa;
  x = inner_product_ab * inv_inner_product_aa;
  return inv_inner_product_aa;
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
