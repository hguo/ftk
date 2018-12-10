#ifndef _FTK_INVERSE_LINEAR_INTERPOLATION_SOLVER_HH
#define _FTK_INVERSE_LINEAR_INTERPOLATION_SOLVER_HH

#include <ftk/numeric/linear_solver.hh>

namespace ftk {

template <typename T>
bool inverse_linear_interpolation1(T v[2], T lambda[1])
{
  return 0;
}

template <typename T>
bool inverse_linear_interpolation3_2(const T u[3], const T v[3], T lambda[3])
{
  const T A[2][2] = {
    {u[0] - u[2], u[1] - u[2]}, 
    {v[0] - v[2], v[1] - v[2]}};
  const T b[2] = {-u[2], -v[2]};

  linear_solver2(A, b, lambda);
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

}

#endif
