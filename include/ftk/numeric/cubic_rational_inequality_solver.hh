#ifndef _FTK_CUBIC_RATIONAL_INEQUALITY_SOLVER_HH
#define _FTK_CUBIC_RATIONAL_INEQUALITY_SOLVER_HH

#include <ftk/numeric/cubic_inequality_solver.hh>
#include <ftk/numeric/rational_function.hh>
#include <map>
#include <set>

namespace ftk {

// returns all intervals that P(x)/Q(x)>=0
template <typename T>
std::tuple<disjoint_intervals<long long>, std::map<long long, T> >
solve_cubic_rational_inequality_quantized(
    const T P[4], const T Q[4], const long long factor = 1000000000L, const T epsilon=std::numeric_limits<T>::epsilon())
{
  disjoint_intervals<long long> I;
  std::map<long long, T> quantized_roots;
  
  if (polynomial_equals_to_zero(P, 3)) { // if the numerator is zero, there is no need to compute the intervals.
    I.set_to_complete();
    return std::make_tuple(I, quantized_roots);
  } else if (polynomial_equals_to_zero(Q, 3)) { // if the denominator is zero, there is no need to compute the intervals as well.
    I.set_to_empty();
    return std::make_tuple(I, quantized_roots);
  }

  // real roots of P and Q
  const auto real_roots_p = solve_cubic_real_multiplicity(P, epsilon), 
             real_roots_q = solve_cubic_real_multiplicity(Q, epsilon);

  // quantized roots of P and Q;  critical values
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

#if 0
  std::cerr << "roots of Q: ";
  for (const auto kv : real_roots_q) 
    std::cerr << "(" << kv.first << "," << kv.second << ")";
  std::cerr << std::endl;
  
  std::cerr << "roots of P: ";
  for (const auto kv : real_roots_p) 
    std::cerr << "(" << kv.first << "," << kv.second << ")";
  std::cerr << std::endl;
#endif

  // singularities
  std::set<long long> singularities;
  for (auto &kv : quantized_roots_q) {
    if (quantized_roots_p.find(kv.first) != quantized_roots_p.end()) { // remove "canceled" roots of Q
      while (quantized_roots_p[kv.first] > 0 && quantized_roots_q[kv.first] > 0) {
        quantized_roots_p[kv.first] --;
        quantized_roots_q[kv.first] --;
        // fprintf(stderr, "reducing multiple roots %lld\n", kv.first);
      }
    }
    if (kv.second > 0) {
      // fprintf(stderr, "adding singularity %lld\n", kv.first);
      singularities.insert(kv.first);
    }
  }

#if 0
  fprintf(stderr, "singularities:\n");
  for (auto s : singularities)
    fprintf(stderr, "**%lld\n", s);
#endif

  // adding quantized roots of P and Q to critical values
  std::set<long long> critical_values;
  for (const auto kv : quantized_roots_p) {
    if (kv.second > 0)
      critical_values.insert(kv.first);
  }
  for (const auto kv : quantized_roots_q) {
    if (kv.second > 0)
      critical_values.insert(kv.first);
  }

  // adding infinities to critical values
  critical_values.insert(basic_interval<long long>::lower_inf());
  critical_values.insert(basic_interval<long long>::upper_inf());

  // sort critical values
  std::vector<long long> sorted_critical_values;
  for (auto v : critical_values)
    sorted_critical_values.push_back(v);
  // std::sort(critical_values.begin(), critical_values.end());
 
#if 0
  fprintf(stderr, "sorted critical values:\n");
  for (int i = 0; i < sorted_critical_values.size(); i ++)
    fprintf(stderr, "##%lld\n", sorted_critical_values[i]);
#endif
  
  // adding subintervals
  for (int i = 0; i < sorted_critical_values.size() - 1; i ++) { // the size is at least 2
    const auto v0 = sorted_critical_values[i], v1 = sorted_critical_values[i+1];
    basic_interval<long long> ii(v0, v1);
    
    if (singularities.find(v0) != singularities.end()) ii.set_lower_open();
    if (singularities.find(v1) != singularities.end()) ii.set_upper_open();
        
    const T x = T(ii.sample()) / factor; // FIXME: avoid sampling on Q's roots
    const T y = evaluate_rational(P, Q, 3, x);
    // std::cerr << i << "checking interval: " << ii << std::endl;
    // std::cerr << "sample: " << x << ", value: " << y << std::endl;
    if (y >= T(0)) {
      I.join(ii);
      // std::cerr << "joining interval " << ii << " I=" << I << std::endl;
    }
  }

  return std::make_tuple(I, quantized_roots);
}

template <typename T>
disjoint_intervals<T> solve_cubic_rational_inequality(
    const T P[4], const T Q[4], const long long factor = 1000000000L, const T epsilon=std::numeric_limits<T>::epsilon())
{
  const auto I0 = std::get<0>(solve_cubic_rational_inequality_quantized(P, Q, factor, epsilon)); // TODO
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
