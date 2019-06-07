#ifndef _FTK_INVERSE_LINEAR_INTERPOLATION_SOLVER_HH
#define _FTK_INVERSE_LINEAR_INTERPOLATION_SOLVER_HH

#include <ftk/numeric/linear_solver.hh>
#include <ftk/numeric/cond.hh>
#include <ftk/numeric/print.hh>

namespace ftk {

template <typename T>
bool inverse_lerp_s2v2(const T V[3][2], T mu[3], const T epsilon = std::numeric_limits<T>::epsilon())
{
  const T A[2][2] = {
    {V[0][0] - V[2][0], V[1][0] - V[2][0]},
    {V[0][1] - V[2][1], V[1][1] - V[2][1]}
  };
  const T b[2] = {-V[2][0], -V[2][1]};

  solve_linear2x2(A, b, mu);
  mu[2] = T(1) - mu[0] - mu[1];

  return mu[0] >= -epsilon && mu[0] <= T(1) + epsilon &&
    mu[1] >= -epsilon && mu[1] <= T(1) + epsilon &&
    mu[2] >= -epsilon && mu[2] <= T(1) + epsilon;
}

template <typename T>
T cond_inverse_lerp_s2v2(const T V[3][2])
{
  const T A[2][2] = {
    {V[0][0] - V[2][0], V[1][0] - V[2][0]},
    {V[0][1] - V[2][1], V[1][1] - V[2][1]}
  };
  return cond_real2x2(A);
}

template <typename T>
bool inverse_lerp_s3v3(const T V[4][3], T lambda[4])
{
  const T A[3][3] = { // linear transformation
    {V[0][0] - V[3][0], V[1][0] - V[3][0], V[2][0] - V[3][0]}, 
    {V[0][1] - V[3][1], V[1][1] - V[3][1], V[2][1] - V[3][1]}, 
    {V[0][2] - V[3][2], V[1][2] - V[3][2], V[2][2] - V[3][2]}
  };
  const T b[3] = {-V[3][0], -V[3][1], -V[3][2]};

  solve_linear3x3(A, b, lambda);
  lambda[3] = T(1) - lambda[0] - lambda[1] - lambda[2];
  
  return lambda[0] >= T(0) && lambda[0] < T(1) && 
         lambda[1] >= T(0) && lambda[1] < T(1) && 
         lambda[2] >= T(0) && lambda[2] < T(1) && 
         lambda[3] >= T(0) && lambda[3] < T(1);
}

}

#endif
