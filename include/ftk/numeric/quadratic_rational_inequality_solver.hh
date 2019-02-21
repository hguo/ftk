#ifndef _FTK_QUADRATIC_RATIONAL_INEQUALITY_SOLVER_HH
#define _FTK_QUADRATIC_RATIONAL_INEQUALITY_SOLVER_HH

#include <ftk/numeric/quadratic_solver.hh>
#include <ftk/numeric/rational.hh>
#include <ftk/numeric/disjoint_intervals.hh>

namespace ftk {

template <typename T>
std::tuple<disjoint_intervals<long long>, std::map<long long, T>>
solve_quadratic_rational_inequality_quantized(
    const T P[3], const T Q[3], const long long factor = 1000000000L, 
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  disjoint_intervals<long long> I;
  std::map<long long, T> quantized_roots;

  if (polynomial_equals_to_zero(P, 2)) {
    I.set_to_complete();
    return std::make_tuple(I, quantized_roots);
  } else if (polynomial_equals_to_zero(Q, 3)) {
    I.set_to_empty();
    return std::make_tuple(I, quantized_roots);
  }

  const auto real_roots_p = solve_quadratic_real_multiplicity(P, epsilon), 
             real_roots_q = solve_quadratic_real_multiplicity(Q, epsilon);

  std::map<long long, int> quantized_roots_p, quantized_roots_q;

  for (const auto kv : real_roots_p) {
    const long long v = kv.first * T(factor);
    quantized_roots_p[v] = kv.second;
    quantized_roots[v] = kv.first;
  }
  for (const auto kv : real_roots_q) {
    const long long v = kv.first * T(factor);
    quantized_roots_q[v] = kv.second;
    quantized_roots[v] = kv.first;
  }

  // singularities
  std::set<long long> singularities;
  for (auto &kv : quantized_roots_q) {
    if (quantized_roots_p.find(kv.first) != quantized_roots_p.end()) {
      while (quantized_roots_p[kv.first] > 0 && quantized_roots_q[kv.first] > 0) {
        quantized_roots_p[kv.first] --;
        quantized_roots_q[kv.first] --;
      }
    }
    if (kv.second > 0) 
      singularities.insert(kv.first);
  }

  std::set<long long> critical_values;
  for (const auto kv : quantized_roots_p) {
    if (kv.second > 0)
      critical_values.insert(kv.first);
  }
  for (const auto kv : quantized_roots_q) {
    if (kv.second > 0) 
      critical_values.insert(kv.first);
  }

  critical_values.insert(basic_interval<long long>::lower_inf());
  critical_values.insert(basic_interval<long long>::upper_inf());

  std::vector<long long> sorted_critical_values;
  for (auto v : critical_values)
    sorted_critical_values.push_back(v);

  for (int i = 0; i < sorted_critical_values.size() - 1; i ++) {
    const auto v0 = sorted_critical_values[i], v1 = sorted_critical_values[i+1];
    basic_interval<long long> ii(v0, v1);

    if (singularities.find(v0) != singularities.end()) ii.set_lower_open();
    if (singularities.find(v1) != singularities.end()) ii.set_upper_open();

    const T x = ii.sample() / factor;
    const T y = evaluate_rational(P, Q, 2, x);

    if (y >= T(0)) 
      I.join(ii);
  }

  return std::make_tuple(I, quantized_roots);
}

}

#endif
