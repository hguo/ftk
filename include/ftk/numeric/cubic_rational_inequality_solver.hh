#ifndef _FTK_CUBIC_RATIONAL_INEQUALITY_SOLVER_HH
#define _FTK_CUBIC_RATIONAL_INEQUALITY_SOLVER_HH

#include <ftk/numeric/cubic_inequality_solver.hh>

// returns all intervals that P(x)/Q(x)>=0
namespace <typename T>
inline disjoint_interval<T> solve_cubic_rational_inequality(
    const T P[], const T Q[], const T epsilon=1e-9)
{
  T p[3], q[3]; // roots of P and Q, respectively
  const np = solve_cubic_real
}

#endif
