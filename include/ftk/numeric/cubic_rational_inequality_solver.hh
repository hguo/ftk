#ifndef _FTK_CUBIC_RATIONAL_INEQUALITY_SOLVER_HH
#define _FTK_CUBIC_RATIONAL_INEQUALITY_SOLVER_HH

#include <ftk/numeric/cubic_inequality_solver.hh>
#include <ftk/numeric/rational.hh>

// returns all intervals that P(x)/Q(x)>=0
namespace ftk {
  
template <typename T>
disjoint_intervals<long long> solve_cubic_rational_inequality_quantized(
    const T P[3], const T Q[3], const long long factor = 1000000000L)
{
  const T epsilon = T(1) / T(factor);
  disjoint_intervals<long long> I;

  T p[3] = {0}, q[3] = {0}; // roots of P and Q, respectively
  const int np = solve_cubic_real(P, p, epsilon);
  const int nq = solve_cubic_real(Q, q, epsilon);
  const int n_roots = np + nq;

  std::set<std::pair<long long, bool> > roots; // the boolean value indicates whether the root is from Q
  
  for (int i = 0; i < np; i ++)
    roots.insert(std::make_pair(p[i] * T(factor), false));
  for (int i = 0; i < nq; i ++)
    roots.insert(std::make_pair(q[i] * T(factor), true));

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
        basic_interval<long long> ii(basic_interval<long long>::lower_inf(), it->first);
        if (it->second)
          ii.set_upper_open();
        subintervals.push_back(ii);
        // fprintf(stderr, "First interval, root=%lld\n", *it);
        // std::cerr << subintervals[0] << std::endl;
      }
      else if (i == roots.size() - 1) { // last root
        basic_interval<long long> ii(it->first, basic_interval<long long>::upper_inf());
        if (it->second)
          ii.set_lower_open();
        subintervals.push_back(ii);
        break;
      }
      else {
        auto it1 = std::next(it);
        basic_interval<long long> ii(it->first, it1->first);
        if (it->second)
          ii.set_lower_open();
        if (it1->second)
          ii.set_upper_open();
        subintervals.push_back(ii);
        // subintervals.push_back(basic_interval<long long>(it->first, (std::next(it)->first)));
      }
      i ++;
    }
    // for (auto i : subintervals)
    //   std::cerr << i << std::endl;

    for (auto i : subintervals) {
      const T x = i.sample() * epsilon;
      const T y = rational_evaluate(P, 3, Q, 3, x);
      // std::cerr << "Interval: " << i << std::endl;
      // std::cerr << "Sample: " << x << ", Value: " << y << std::endl;
      if (y >= T(0)) I.join(i);
    }
  }

  // std::cerr << I << std::endl;
  return I;
}

template <typename T>
disjoint_intervals<T> solve_cubic_rational_inequality(
    const T P[], const T Q[], const long long factor = 1000000000L)
{
  const auto I0 = solve_cubic_rational_inequality_quantized(P, Q, factor);
  return disjoint_intervals<T>(I0, factor);

#if 0
  const T epsilon = T(1) / T(factor);
  disjoint_intervals<T> I;
  for (auto i : I0.subintervals()) { // TODO: consider open intervals
    const auto lb = i.lower_bounded() ? (i.lower() * epsilon) : (basic_interval<T>::lower_inf()), 
               ub = i.upper_bounded() ? (i.upper() * epsilon) : (basic_interval<T>::upper_inf());
    I.join(basic_interval<T>(lb, ub));
  }

  return I;
#endif
}

}

#endif
