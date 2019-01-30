#ifndef _FTK_CUBIC_RATIONAL_INEQUALITY_SOLVER_HH
#define _FTK_CUBIC_RATIONAL_INEQUALITY_SOLVER_HH

#include <ftk/numeric/cubic_inequality_solver.hh>
#include <ftk/numeric/rational.hh>

// returns all intervals that P(x)/Q(x)>=0
namespace ftk {
  
template <typename T>
inline disjoint_intervals<long long> solve_cubic_rational_inequality_quntized(
    const T P[], const T Q[], const T factor = 1000000000L)
{
  const T epsilon = T(1) / T(factor);

  T p[3], q[3]; // roots of P and Q, respectively
  const int np = solve_cubic_real(P, p, epsilon), 
            nq = solve_cubic_real(Q, q, epsilon);
  const int n_roots = np + nq;

  // fprintf(stderr, "np=%d, roots={%f, %f, %f}\n", 
  //     np, p[0], p[1], p[2]);
  // fprintf(stderr, "nq=%d, roots={%f, %f, %f}\n", 
  //     nq, q[0], q[1], q[2]);

  std::set<long long> roots;
  
  for (int i = 0; i < np; i ++)
    roots.insert(p[i] * T(factor));
  for (int i = 0; i < nq; i ++)
    roots.insert(q[i] * T(factor));

  disjoint_intervals<long long> I;
  if (roots.size() == 0) { // neither P or Q has roots
    if (std::abs(Q[0]) < epsilon)
      I.set_to_empty();
    else if (P[0] / Q[0] >= T(0))
      I.set_to_complete();
  } else {
    int i = 0;
    std::vector<basic_interval<long long> > subintervals;

    for (auto it = roots.begin(); it != roots.end(); it ++) {
      if (i == 0) { // first root
        subintervals.push_back(basic_interval<long long>(basic_interval<long long>::lower_inf(), *it));
        // fprintf(stderr, "First interval, root=%lld\n", *it);
        // std::cerr << subintervals[0] << std::endl;
      }
      else if (i == roots.size() - 1) { // last root
        subintervals.push_back(basic_interval<long long>(*it, basic_interval<long long>::upper_inf()));
        break;
      }
      else subintervals.push_back(basic_interval<long long>(*it, *(std::next(it))));
      i ++;
    }

    // fprintf(stderr, "roots:\n");
    // for (auto r : roots) {
    //   fprintf(stderr, "%lld\n", r);
    // }

    for (auto i : subintervals) {
      const T x = i.sample() * epsilon;
      const T y = rational_evaluate(P, 3, Q, 3, x);
      // std::cerr << "Interval: " << i << std::endl;
      // std::cerr << "Sample: " << x << ", Value: " << y << std::endl;
      if (y >= T(0)) I.join(i);
    }
  }

  return I;
}

template <typename T>
inline disjoint_intervals<T> solve_cubic_rational_inequality(
    const T P[], const T Q[], const T factor = 1000000000L)
{
  const auto I0 = solve_cubic_rational_inequality_quntized(P, Q, factor);
  const T epsilon = T(1) / T(factor);

  disjoint_intervals<T> I;
  for (auto i : I0.subintervals()) {
    const auto lb = i.lower_bounded() ? (i.lower() * epsilon) : (basic_interval<T>::lower_inf()), 
               ub = i.upper_bounded() ? (i.upper() * epsilon) : (basic_interval<T>::upper_inf());
    I.join(basic_interval<T>(lb, ub));
  }

  return I;
}

}

#endif
