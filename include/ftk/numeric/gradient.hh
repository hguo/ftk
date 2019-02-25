#ifndef _FTK_GRADIENT_HH
#define _FTK_GRADIENT_HH

#include <ftk/numeric/linear_solver.hh>

namespace ftk {

template <typename T>
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
void gradient_2dsimplex2_2(const T X[3][2], const T f[3][2], T gradf[2][2])
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

}

#endif
