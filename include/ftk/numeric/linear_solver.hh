#ifndef _FTK_LINEAR_SOLVE3_H
#define _FTK_LINEAR_SOLVE3_H

#include <ftk/numeric/det.hh>

namespace ftk {

template <typename T>
inline T linear_solver2(const T A[2][2], const T b[2], T x[2])
{
  const T D  = A[0][0]*A[1][1] - A[0][1]*A[1][0],
          Dx = b[0]   *A[1][1] - A[0][1]*b[1],
          Dy = A[0][0]*b[1]    - b[0]   *A[1][0];

  x[0] = Dx / D;
  x[1] = Dy / D;
  
  return D;
}

}

#endif
