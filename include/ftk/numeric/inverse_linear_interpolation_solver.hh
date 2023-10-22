#ifndef _FTK_INVERSE_LINEAR_INTERPOLATION_SOLVER_HH
#define _FTK_INVERSE_LINEAR_INTERPOLATION_SOLVER_HH

#include <ftk/config.hh>
#include <ftk/numeric/linear_solver.hh>
#include <ftk/numeric/cond.hh>
#include <ftk/numeric/print.hh>

namespace ftk {

template <typename T>
__device__ __host__
inline bool inverse_lerp_s1v1(const T f[2], T mu[2])
{
  mu[0] = f[1] / (f[1] - f[0]);
  mu[1] = T(1) - mu[0];

  return mu[0] >= T(0) && mu[0] < T(1);
}

template <typename T>
__device__ __host__
inline bool inverse_lerp_s1v2(const T V[2][2], T &mu, const T epsilon = std::numeric_limits<T>::epsilon())
{
  const T a[2] = {V[0][0] - V[1][0], V[0][1] - V[1][1]}, 
          b[2] = {-V[1][0], -V[1][1]};

  T res = solve_least_square2x1(a, b, mu);
  return mu >= T(0) && mu < T(1); // && std::abs(res) < 0.01; //epsilon;
}

template <typename T>
__device__ __host__
inline bool inverse_lerp_s2v2(const T V[3][2], T mu[3], 
    T *cond = NULL,
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  const T A[2][2] = {
    {V[0][0] - V[2][0], V[1][0] - V[2][0]},
    {V[0][1] - V[2][1], V[1][1] - V[2][1]}
  };

  if (cond)
    *cond = cond_real2x2(A);

  const T b[2] = {-V[2][0], -V[2][1]};

  solve_linear2x2(A, b, mu);
  mu[2] = T(1) - mu[0] - mu[1];

  return mu[0] >= -epsilon && mu[0] <= T(1) + epsilon &&
    mu[1] >= -epsilon && mu[1] <= T(1) + epsilon &&
    mu[2] >= -epsilon && mu[2] <= T(1) + epsilon;
}

template <typename T>
__device__ __host__
inline bool inverse_lerp_s2v3(const T V[3][3], T mu[3], // will use least squares to figure out the barycentric coordinates
    T *cond = NULL,
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  const T A[3][2] = {
    {V[0][0] - V[2][0], V[1][0] - V[2][0]},
    {V[0][1] - V[2][1], V[1][1] - V[2][1]},
    {V[0][2] - V[2][2], V[1][2] - V[2][2]}
  };
  // print3x2("A", A);

  const T b[3] = {-V[2][0], -V[2][1], -V[2][2]};
  // print3("b", b);

  T cond_ = solve_least_square3x2(A, b, mu);
  // print3("mu0", mu);

  if (cond)
    *cond = cond_;
  
  mu[2] = T(1) - mu[0] - mu[1];

  return mu[0] >= -epsilon && mu[0] <= T(1) + epsilon &&
    mu[1] >= -epsilon && mu[1] <= T(1) + epsilon &&
    mu[2] >= -epsilon && mu[2] <= T(1) + epsilon;
}

template <typename T>
__device__ __host__
inline bool inverse_lerp_s2v3_0(const T V[3][3], const T V0[3], // target value
    T mu[3], // barycentric
    T *cond = NULL,
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  const T W[3][3] = {
    {V[0][0] - V0[0], V[0][1] - V0[1], V[0][2] - V0[2]},
    {V[1][0] - V0[0], V[1][1] - V0[1], V[1][2] - V0[2]},
    {V[2][0] - V0[0], V[2][1] - V0[1], V[2][2] - V0[2]}
  };
  // print3x3("W", W);
  return inverse_lerp_s2v3(W, mu, cond, epsilon);
}

template <typename T>
__device__ __host__
T cond_inverse_lerp_s2v2(const T V[3][2])
{
  const T A[2][2] = {
    {V[0][0] - V[2][0], V[1][0] - V[2][0]},
    {V[0][1] - V[2][1], V[1][1] - V[2][1]}
  };
  return cond_real2x2(A);
}

template <typename L>
__device__ __host__
inline bool integer_inverse_lerp_s1v2(const L V[2][2])
{
  const L A[2] = {V[0][0] - V[1][0], V[0][1] - V[1][1]},
          B[2] = {-V[1][0], -V[1][1]};

  int ns0 = number_roots_linear1(A[0], B[0]),
      ns1 = number_roots_linear1(A[1], B[1]);

  auto check_range = [](L a, L b) { // check if (b/a) \in (0, 1)
    if (sign(a) != sign(b)) return false;
    a = std::abs(a);
    b = std::abs(b);
    return (b > 0) && (b < a);
  };

  if (ns0 == LINEAR_EQUATION_NO_SOLUTION || ns1 == LINEAR_EQUATION_NO_SOLUTION) // no solution at all
    return false;
  else if (ns0 == LINEAR_EQUATION_UNIQUE_SOLUTION && ns1 == LINEAR_EQUATION_UNIQUE_SOLUTION) // both equations have identical solution
    if (A[0] * B[1] == A[1] * B[0]) 
      return check_range(A[0], B[0]);
    else 
      return false;
  else if (ns0 == LINEAR_EQUATION_INFINITE_SOLUTIONS && ns1 == LINEAR_EQUATION_UNIQUE_SOLUTION) 
    return check_range(A[1], B[1]);
  else if (ns0 == LINEAR_EQUATION_UNIQUE_SOLUTION && ns1 == LINEAR_EQUATION_INFINITE_SOLUTIONS) 
    return check_range(A[0], B[0]);
  else return false; // both equations have infinite solutions
}

template <typename T>
__device__ __host__
inline bool inverse_lerp_s3v3(const T V[4][3], T lambda[4],
    T *cond = NULL,
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  const T A[3][3] = { // linear transformation
    {V[0][0] - V[3][0], V[1][0] - V[3][0], V[2][0] - V[3][0]}, 
    {V[0][1] - V[3][1], V[1][1] - V[3][1], V[2][1] - V[3][1]}, 
    {V[0][2] - V[3][2], V[1][2] - V[3][2], V[2][2] - V[3][2]}
  };

  if (cond) 
    *cond = cond_real3x3(A);

  const T b[3] = {-V[3][0], -V[3][1], -V[3][2]};

  solve_linear3x3(A, b, lambda);
  lambda[3] = T(1) - lambda[0] - lambda[1] - lambda[2];
  
  return lambda[0] >= -epsilon && lambda[0] < T(1) + epsilon && 
         lambda[1] >= -epsilon && lambda[1] < T(1) + epsilon && 
         lambda[2] >= -epsilon && lambda[2] < T(1) + epsilon && 
         lambda[3] >= -epsilon && lambda[3] < T(1) + epsilon;
}

}

#endif
