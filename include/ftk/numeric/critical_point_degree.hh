#ifndef _FTK_CRITICAL_POINT_DEGREE_HH
#define _FTK_CRITICAL_POINT_DEGREE_HH

#include <ftk/config.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/sign.hh>

namespace ftk {

template <typename T>
int critical_point_degree_simplex2(const T V[3][2], const T X[3][3])
{
  const T m0 = vector_2norm_2(V[0]), 
          m1 = vector_2norm_2(V[1]),
          m2 = vector_2norm_2(V[2]);

  const T cross01 = cross_product2(V[0], V[1]) / (m0 * m1), 
          cross12 = cross_product2(V[1], V[2]) / (m1 * m2),
          cross20 = cross_product2(V[2], V[0]) / (m2 * m0);

  int cnt_neg = 0;
  if (cross01 < 0) cnt_neg ++;
  if (cross12 < 0) cnt_neg ++;
  if (cross20 < 0) cnt_neg ++;

  int orientation = sign((X[1][1] - X[0][1]) * (X[2][0] - X[1][0]) - (X[2][1] - X[1][1]) * (X[1][0] - X[0][0]));
  int deg = cnt_neg >= 2 ? -1 : 1;

  return deg * orientation;
}

}

#endif
