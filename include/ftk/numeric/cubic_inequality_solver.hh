#ifndef _FTK_CUBIC_INEQUALITY_SOLVER_HH
#define _FTK_CUBIC_INEQUALITY_SOLVER_HH

#include <ftk/numeric/disjoint_intervals.hh>
#include <ftk/numeric/cubic_solver.hh>
#include <ftk/numeric/polynomial.hh>
#include <vector>

namespace ftk {

// returns all intervals that P(x)>=0
template <typename T>
inline disjoint_intervals<T> solve_cubic_inequality_real(
    const T P[4], const T epsilon=std::numeric_limits<T>::epsilon())
{
  T x[3]; // roots
  const int n_roots = solve_cubic_real(P, x, epsilon);

  if (n_roots == 0) {
    disjoint_intervals<T> I;
    if (P[0] >= T(0)) I.set_to_complete();
    else I.set_to_empty();
    return I;
  } else {
    std::vector<basic_interval<T> > subintervals;
    std::sort(x, x+n_roots);
    if (n_roots == 3) {
      subintervals.push_back(basic_interval<T>(basic_interval<T>::lower_inf(), x[0]));
      subintervals.push_back(basic_interval<T>(x[0], x[1]));
      subintervals.push_back(basic_interval<T>(x[1], x[2]));
      subintervals.push_back(basic_interval<T>(x[2], basic_interval<T>::upper_inf()));
    } else if (n_roots == 2) {
      subintervals.push_back(basic_interval<T>(basic_interval<T>::lower_inf(), x[0]));
      subintervals.push_back(basic_interval<T>(x[0], x[1]));
      subintervals.push_back(basic_interval<T>(x[1], basic_interval<T>::upper_inf()));
    } else if (n_roots == 1) {
      subintervals.push_back(basic_interval<T>(basic_interval<T>::lower_inf(), x[0]));
      subintervals.push_back(basic_interval<T>(x[0], basic_interval<T>::upper_inf()));
    }

    disjoint_intervals<T> I;
    for (auto i : subintervals) {
      const T x1 = i.sample();
      const T y = polynomial_evaluate(P, 3, x1);
      // fprintf(stderr, "x=%f, y=%f\n", x1, y);
      if (y >= T(0)) I.join(i);
    }

    return I;
  }
}

}

#endif
