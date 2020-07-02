#ifndef _FTK_SUBMATRIX_HH
#define _FTK_SUBMATRIX_HH

#include <ftk/ftk_config.hh>

namespace ftk {

template <typename T>
void submatrix_3x3_2x2(const T A[3][3], T B[2][2])
{
  B[0][0] = A[0][0];
  B[0][1] = A[0][1];
  B[1][0] = A[1][0];
  B[1][1] = A[1][1];
}

}
