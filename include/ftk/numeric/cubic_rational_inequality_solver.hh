#ifndef _FTK_CUBIC_RATIONAL_INEQUALITY_SOLVER_HH
#define _FTK_CUBIC_RATIONAL_INEQUALITY_SOLVER_HH

#include <ftk/numeric/cubic_inequality_solver.hh>
#include <ftk/numeric/rational.hh>

namespace ftk {

// returns all intervals that P(x)/Q(x)>=0
template <typename T>
disjoint_intervals<long long> solve_cubic_rational_inequality_quantized(
    const T P[4], const T Q[4], const long long factor = 1000000000L)
{
  const T epsilon = T(1) / T(factor);
  disjoint_intervals<long long> I;
  
  if (polynomial_equals_to_zero(P, 3)) { // if the numerator is zero, there is no need to compute the intervals.
    I.set_to_complete();
    return I;
  } else if (polynomial_equals_to_zero(Q, 3)) { // if the denominator is zero, there is no need to compute the intervals as well.
    I.set_to_empty();
    return I;
  }

  T p[3] = {0}, q[3] = {0}; // roots of P and Q, respectively
  const int np = solve_cubic_real(P, p, epsilon);
  const int nq = solve_cubic_real(Q, q, epsilon);

  // fprintf(stderr, "P={%f, %f, %f, %f}, np=%d, roots={%f, %f, %f}\n", 
  //     P[0], P[1], P[2], P[3], np, p[0], p[1], p[2]);
  // fprintf(stderr, "Q={%f, %f, %f, %f}, nq=%d, roots={%f, %f, %f}\n", 
  //     Q[0], Q[1], Q[2], Q[3], nq, q[0], q[1], q[2]);

  // critical values
  std::vector<long long> critical_values;
  for (int i = 0; i < np; i ++)
    critical_values.push_back(p[i] * T(factor));
  for (int i = 0; i < np; i ++)
    critical_values.push_back(q[i] * T(factor));

  // adding infinities to critical values
  critical_values.push_back(basic_interval<long long>::lower_inf());
  critical_values.push_back(basic_interval<long long>::upper_inf());
  
  // sort critical values
  std::sort(critical_values.begin(), critical_values.end());

  // singularities
  std::set<long long> singularities;
  for (int i = 0; i < nq; i ++) {
    if (!isfinite(evaluate_rational(P, Q, 3, q[i])))
      singularities.insert(q[i] * T(factor));
  }
  
  // adding subintervals
  for (int i = 0; i < critical_values.size() - 1; i ++) { // the size is at least 2
    const auto v0 = critical_values[i], v1 = critical_values[i+1];
    basic_interval<long long> ii(v0, v1);
    
    if (singularities.find(v0) != singularities.end()) ii.set_lower_open();
    if (singularities.find(v1) != singularities.end()) ii.set_upper_open();
        
    // std::cerr << i << "checking interval: " << ii << std::endl;
    
    const T x = ii.sample() * epsilon; // FIXME: avoid sampling on Q's roots
    const T y = evaluate_rational(P, Q, 3, x);
    // std::cerr << "sample: " << x << ", value: " << y << std::endl;
    if (y >= T(0)) I.join(ii);
  }
  
  return I;
}

template <typename T>
disjoint_intervals<T> solve_cubic_rational_inequality(
    const T P[4], const T Q[4], const long long factor = 1000000000L)
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
