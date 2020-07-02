#ifndef _FTK_CHARACTERISTIC_POLYNOMIAL_HH
#define _FTK_CHARACTERISTIC_POLYNOMIAL_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/trace.hh>
#include <ftk/numeric/det.hh>

namespace ftk {

template <typename T>
__host__ __device__
void characteristic_polynomial_2x2(const T A[2][2], T P[3])
{
  P[2] = T(1);
  P[1] = -trace2(A);
  P[0] = det2(A);
}

template <typename T>
__host__ __device__
void characteristic_polynomial_2x2(const T A[2][2], const T B[2][2], T P[3])
{
  P[2] = det2(B);
  P[1] = -(A[0][0]*B[1][1] - A[1][0]*B[0][1] + B[0][0]*A[1][1] - B[1][0]*A[0][1]);
  P[0] = det2(A);
}

template <typename T>
__host__ __device__
void characteristic_polynomial_2x2(T a00, T a01, T a10, T a11, T b00, T b01, T b10, T b11, T P[3])
{
  P[2] = b00 * b11 - b10 * b01;
  P[1] = -(a00 * b11 - a10 * b01 + b00 * a11 - b10 * a01);
  P[0] = a00 * a11 - a10 * a01;
}

// used for eigensolver
template <typename T>
__host__ __device__
void characteristic_polynomial_3x3(const T A[3][3], T P[4])
{
  P[3] = 1;
  P[2] = -trace3(A);
  P[1] = A[1][1]*A[2][2] + A[0][0]*A[2][2] + A[0][0]*A[1][1] 
    - A[0][1]*A[1][0] - A[1][2]*A[2][1] - A[0][2]*A[2][0];
  P[0] = -det3(A);
}

// used for generalized eigensolver
template <typename T>
__host__ __device__
void characteristic_polynomial_3x3(const T A[3][3], const T B[3][3], T P[4]) // tested.
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
  P[3] = -det3(B);
}

}

#endif
