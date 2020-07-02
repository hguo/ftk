#ifndef _FTK_GRADIENT_HH
#define _FTK_GRADIENT_HH

#include <ftk/numeric/linear_solver.hh>

namespace ftk {

template <typename T>
__device__ __host__
void gradient_2dsimplex2(const T X[3][2], const T f[3], T gradf[2])
{
  const T A[3][2] = {
    {X[1][0] - X[0][0], X[1][1] - X[0][1]}, 
    {X[2][0] - X[1][0], X[2][1] - X[1][1]},
    {X[0][0] - X[2][0], X[0][1] - X[2][1]}
  };
  const T b[3] = {f[1] - f[0], f[2] - f[1], f[0] - f[2]};

  solve_least_square3x2(A, b, gradf);
}

template <typename T>
__device__ __host__
void jacobian_2dsimplex2(const T X[3][2], const T f[3][2], T gradf[2][2])
// gradient_2dsimplex2_2
{
  const T A[3][2] = {
    {X[1][0] - X[0][0], X[1][1] - X[0][1]}, 
    {X[2][0] - X[1][0], X[2][1] - X[1][1]},
    {X[0][0] - X[2][0], X[0][1] - X[2][1]}
  };
  const T B[3][2] = {
    {f[1][0] - f[0][0], f[1][1] - f[0][1]},
    {f[2][0] - f[1][0], f[2][1] - f[1][1]},
    {f[0][0] - f[2][0], f[0][1] - f[2][1]}
  };
  
  solve_least_square3x2_2(A, B, gradf);
}

template <typename T>
__device__ __host__
void jacobian_3dsimplex(const T X[4][3], const T V[4][3], T J[3][3])
// gradient_2dsimplex2_2
{
  const T A[4][3] = {
    {X[1][0] - X[0][0], X[1][1] - X[0][1], X[1][2] - X[0][2]},
    {X[2][0] - X[1][0], X[2][1] - X[1][1], X[2][2] - X[1][2]},
    {X[3][0] - X[2][0], X[3][1] - X[2][1], X[3][2] - X[2][2]},
    {X[0][0] - X[3][0], X[0][1] - X[3][1], X[0][2] - X[3][2]}
  };

  const T B[4][3] = {
    {V[1][0] - V[0][0], V[1][1] - V[0][1], V[1][2] - V[0][2]},
    {V[2][0] - V[1][0], V[2][1] - V[1][1], V[2][2] - V[1][2]},
    {V[3][0] - V[2][0], V[3][1] - V[2][1], V[3][2] - V[2][2]},
    {V[0][0] - V[3][0], V[0][1] - V[3][1], V[0][2] - V[3][2]}
  };
  
  solve_least_square4x3_3(A, B, J);
}

}

#endif
