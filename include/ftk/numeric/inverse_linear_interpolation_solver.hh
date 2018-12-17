#ifndef _FTK_INVERSE_LINEAR_INTERPOLATION_SOLVER_HH
#define _FTK_INVERSE_LINEAR_INTERPOLATION_SOLVER_HH

#include <ftk/numeric/linear_solver.hh>
#include <ftk/numeric/print.hh>

namespace ftk {

template <typename T>
bool inverse_linear_interpolation1(T v[2], T lambda[1])
{
  return 0;
}

template <typename T>
bool inverse_linear_interpolation3_2(const T u[3], const T v[3], T lambda[3], T &cond)
{
  const T A[2][2] = {
    {u[0] - u[2], u[1] - u[2]}, 
    {v[0] - v[2], v[1] - v[2]}};
  const T b[2] = {-u[2], -v[2]};

  // linear_solver2(A, b, lambda);
  cond = linear_solver2_cond(A, b, lambda);
  lambda[2] = T(1) - lambda[0] - lambda[1];

  return lambda[0] >= T(0) && lambda[0] < T(1) && 
         lambda[1] >= T(0) && lambda[1] < T(1) && 
         lambda[2] >= T(0) && lambda[2] < T(1);
}

template <typename T>
bool inverse_linear_interpolation3(T v[4][3], T lambda[3])
{
  return 0;
}

template <typename T>
bool inverse_linear_interpolation4_3(const T V[4][3], T lambda[4])
{
  const T A[3][3] = { // linear transformation
    {V[0][0] - V[3][0], V[1][0] - V[3][0], V[2][0] - V[3][0]}, 
    {V[0][1] - V[3][1], V[1][1] - V[3][1], V[2][1] - V[3][1]}, 
    {V[0][2] - V[3][2], V[1][2] - V[3][2], V[2][2] - V[3][2]}
  };
  const T b[3] = {-V[3][0], -V[3][1], -V[3][2]};

  linear_solver3(A, b, lambda);
  lambda[3] = T(1) - lambda[0] - lambda[1] - lambda[2];
  
  return lambda[0] >= T(0) && lambda[0] < T(1) && 
         lambda[1] >= T(0) && lambda[1] < T(1) && 
         lambda[2] >= T(0) && lambda[2] < T(1) && 
         lambda[3] >= T(0) && lambda[3] < T(1);
}

}

#endif
