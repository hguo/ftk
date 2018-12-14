#ifndef _FTK_CHARACTERISTIC_POLYNOMIAL_HH
#define _FTK_CHARACTERISTIC_POLYNOMIAL_HH

#include <ftk/numeric/trace.hh>
#include <ftk/numeric/det.hh>

namespace ftk {

// used for eigensolver
template <typename T>
void characteristic_polynomial_3x3(const T A[3][3], T P[4])
{
  P[3] = 1;
  P[2] = -trace3(A);
  P[1] = A[1][1]*A[2][2] + A[0][0]*A[2][2] + A[0][0]*A[1][1] 
    - A[0][1]*A[1][0] - A[1][2]*A[2][1] - A[0][2]*A[2][0];
  P[0] = -det3(A);
}

// used for generalized eigensolver
void characteristic_polynomial_3x3(const T A[3][3], const T B[3][3], T P[4])
{
  P[0] = det3(A);
  P[1] = A[1][2] * A[2][1] * B[0][0] - A[1][1] * A[2][2] * B[0][0] - A[1][2] * A[2][0] * B[0][1]  
        +A[1][0] * A[2][2] * B[0][1] + A[1][1] * A[2][0] * B[0][2] - A[1][0] * A[2][1] * B[0][2]
        -A[0][2] * A[2][1] * B[1][0] + A[0][1] * A[2][2] * B[1][0] + A[0][2] * A[2][0] * B[1][1]
        -A[0][0] * A[2][2] * B[1][1] - A[0][1] * A[2][0] * B[1][2] + A[0][0] * A[2][1] * B[1][2]
        +A[0][2] * A[1][1] * B[2][0] - A[0][1] * A[1][2] * B[2][0] - A[0][2] * A[1][0] * B[2][1]
        +A[0][0] * A[1][2] * B[2][1] + A[0][1] * A[1][0] * B[2][2] - A[0][0] * A[1][1] * B[2][2];
  P[2] =-A[2][2] * B[0][1] * B[1][0] + A[2][1] * B[0][2] * B[1][0] + A[2][2] * B[0][0] * B[1][1]
        -A[2][0] * B[0][2] * B[1][1] - A[2][1] * B[0][0] * B[1][2] + A[2][0] * B[0][1] * B[1][2]
        +A[1][2] * B[0][1] * B[2][0] - A[1][1] * B[0][2] * B[2][0] - A[0][2] * B[1][1] * B[2][0]
        +A[0][1] * B[1][2] * B[2][0] - A[1][2] * B[0][0] * B[2][1] + A[1][0] * B[0][2] * B[2][1]
        +A[0][2] * B[1][0] * B[2][1] - A[0][0] * B[1][2] * B[2][1] + A[1][1] * B[0][0] * B[2][2] 
        -A[1][0] * B[0][1] * B[2][2] - A[0][1] * B[1][0] * B[2][2] + A[0][0] * B[1][1] * B[2][2];
  P[3] = det3(B);
}

}

#endif
