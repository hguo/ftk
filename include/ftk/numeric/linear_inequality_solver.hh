#ifndef _FTK_LINEAR_INEQUALITY_SOLVER_HH
#define _FTK_LINEAR_INEQUALITY_SOLVER_HH

#include <ftk/numeric/disjoint_intervals.hh>

namespace ftk {

template <typename T>
inline basic_interval<T> solve_linear_inequality(T a, T b) // ax + b >= 0
{
  basic_interval<T> I;
  if (std::abs(a) < std::numeric_limits<T>::epsilon()) {
    if (b >= T(0)) I.set_to_complete();
    else I.set_to_empty();
  } else if (a > T(0)) {
    I.set_lower(-b / a);
    I.set_lower_closed();
    I.set_upper_inf();
  } else {
    I.set_upper(-b / a);
    I.set_upper_closed();
    I.set_lower_inf();
  }
  return I;
}

}

#endif
