#ifndef _FTK_PARALLEL_VECTOR_SOLVER2_H
#define _FTK_PARALLEL_VECTOR_SOLVER2_H

#include <ftk/ftk_config.hh>
#include <ftk/numeric/matrix_inverse.hh>
#include <ftk/numeric/matrix_addition.hh>
#include <ftk/numeric/eigen_solver2.hh>
#include <ftk/numeric/eigen_solver3.hh>
#include <ftk/numeric/linear_solver.hh>
#include <ftk/numeric/linear_inequality_solver.hh>
#include <ftk/numeric/quadratic_rational_inequality_solver.hh>
#include <ftk/numeric/cubic_rational_inequality_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/transpose.hh>
#include <ftk/numeric/characteristic_polynomial.hh>
#include <ftk/numeric/polynomial.hh>
#include <ftk/numeric/matrix_norm.hh>
#include <ftk/numeric/isnan.hh>
#include <ftk/numeric/print.hh>

namespace ftk {

template <typename T>
inline bool verify_pv_s0v2(const T v[2], const T w[2], const T epsilon = std::numeric_limits<T>::epsilon())
{
  const double c = ftk::cross_product2(v, w); // cross product
  return std::abs(c) <= epsilon;
}

template <typename T>
inline bool verify_pv_s1v2(const T V[2][2], const T W[2][2], const T mu[2], 
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  double v[2], w[2];
  ftk::lerp_s1v2(V, mu, v);
  ftk::lerp_s1v2(W, mu, w);
  return verify_pv_s0v2(v, w, epsilon);
}

template <typename T>
inline bool verify_pv_s2v2(const T V[3][2], const T W[3][2], const T mu[3],
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  double v[2], w[2];
  ftk::lerp_s2v2(V, mu, v);
  ftk::lerp_s2v2(W, mu, w);
  return verify_pv_s0v2(v, w, epsilon);
}

template <typename T>
inline bool verify_pv_s2v2(const T V[3][2], const T w[2], const T mu[3], 
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  double v[2];
  ftk::lerp_s2v2(V, mu, v);
  return verify_pv_s0v2(v, w, epsilon);
}

template <typename T>
inline int solve_pv_s1v2(const T VV[2][2], const T WW[2][2], 
    T lambda[2], T mu[2][2], const T epsilon = std::numeric_limits<T>::epsilon())
{
  T V[2][2], W[2][2];
  transpose2x2(VV, V);
  transpose2x2(WW, W);

  T P[3];
  characteristic_polynomial_2x2(V, W, P);

  T l[2] = {0};
  const int n_roots = solve_quadratic_real(P, l, epsilon);
  int n_solutions = 0;

  for (int i = 0; i < n_roots; i ++) {
    if (std::abs(l[i]) <= epsilon) continue;

    const T a[2] = {
      (V[0][0] - V[0][1]) - l[i] * (W[0][0] - W[0][1]), 
      (V[1][0] - V[1][1]) - l[i] * (W[1][0] - W[1][1])
    };
    const T b[2] = {
      -(V[0][1] - l[i] * W[0][1]),
      -(V[1][1] - l[i] * W[1][1])
    };

    T nu[2];
    const auto cond = solve_least_square2x1(a, b, nu[0]);
    nu[1] = T(1) - nu[0];

    if (nu[0] >= -epsilon && nu[0] <= 1+epsilon && 
        nu[1] >= -epsilon && nu[1] <= 1+epsilon) 
    {
      double v[2], w[2];
      ftk::lerp_s1v2(VV, nu, v);
      ftk::lerp_s1v2(WW, nu, w);
      const double c = ftk::cross_product2(v, w);
      if (c > 1e-2) {
        fprintf(stderr, "rejecting: nu={%f, %f}, v={%f, %f}, w={%f, %f}, c=%f\n", 
            nu[0], nu[1], v[0], v[1], w[0], w[1], c);
        print2x2("V", VV);
        print2x2("W", WW);
        continue;
      }

      const int j = n_solutions ++;
      lambda[j] = l[i];
      mu[j][0] = nu[0];
      mu[j][1] = nu[1];
    }
  }

  return n_solutions;
}

template <typename T>
inline void characteristic_polynomials_pv_s2v2(const T V[3][2], const T w[2], T P[3][2]) // the vector field w is constant
{
  const T A[2][2] = { // linear transformation
    {V[0][0] - V[2][0], V[1][0] - V[2][0]}, 
    {V[0][1] - V[2][1], V[1][1] - V[2][1]}
  };
  T invA[2][2];
  T det = ftk::matrix_inverse2x2(A, invA);

  P[0][0] =-(invA[0][0] * V[2][0] + invA[0][1] * V[2][1]);
  P[0][1] =  invA[0][0] * w[0] + invA[0][1] * w[1];
  P[1][0] =-(invA[1][0] * V[2][0] + invA[1][1] * V[2][1]);
  P[1][1] =  invA[1][0] * w[0] + invA[1][1] * w[1];
  P[2][0] = T(1) - P[0][0] - P[1][0];
  P[2][1] =-(P[0][1] + P[1][1]);
}

template <typename T>
inline void characteristic_polynomials_pv_s2v2(const T V[3][2], const T W[3][2], T Q[3], T P[3][3])
{
  const T A[2][2] = { // linear transformation
    {V[0][0] - V[2][0], V[1][0] - V[2][0]}, 
    {V[0][1] - V[2][1], V[1][1] - V[2][1]}
  };
  const T B[2][2] = {
    {W[0][0] - W[2][0], W[1][0] - W[2][0]}, 
    {W[0][1] - W[2][1], W[1][1] - W[2][1]}
  };
  const T a[2] = {V[2][0], V[2][1]};
  const T b[2] = {W[2][0], W[2][1]};
  
  const T rhs[2][2] = {
    {-a[0], b[0]}, 
    {-a[1], b[1]}
  };

  // coefs for Q
  characteristic_polynomial_2x2(A, B, Q);

  T adj[2][2][2] = { // adjugate matrix
    {{ A[1][1], -B[1][1]}, {-A[0][1],  B[0][1]}}, 
    {{-A[1][0],  B[1][0]}, { A[0][0], -B[0][0]}}
  };

  for (int i = 0; i < 2; i ++) {
    T poly[3] = {0};
    for (int j = 0; j < 3; j ++)
      P[i][j] = 0; // clear results
    for (int j = 0; j < 2; j ++) {
      polynomial_multiplication(adj[i][j], 1, rhs[j], 1, poly);
      polynomial_addition_in_place(P[i], 2, poly, 2);
    }
  }

  P[2][0] = Q[0] - P[0][0] - P[1][0];
  P[2][1] = Q[1] - P[0][1] - P[1][1];
  P[2][2] = Q[2] - P[0][2] - P[1][2];
}

template <typename T>
disjoint_intervals<T> solve_pv_inequalities_s2v2(const T P[3][2], const T w[2])
{
  disjoint_intervals<T> I;
  I.set_to_complete();

  // T P[3][2] = {0};
  // characteristic_polynomials_parallel_vector2_simplex2(V, w, P);
 
  I.intersect( solve_linear_inequality(P[0][1], P[0][0]) );
  I.intersect( solve_linear_inequality(-P[0][1], 1 - P[0][0]) );
  I.intersect( solve_linear_inequality(P[1][1], P[1][0]) );
  I.intersect( solve_linear_inequality(-P[1][1], 1 - P[1][0]) );
  I.intersect( solve_linear_inequality(P[2][1], P[2][0]) );
  I.intersect( solve_linear_inequality(-P[2][1], T(1) - P[2][0]) );

  // std::cerr << I << std::endl;

  return I;
}

template <typename T>
std::tuple<disjoint_intervals<long long>, std::map<long long, T>>
solve_pv_inequalities_quantized_s2v2(
    const T Q[3], const T P[3][3], 
    // const T V[3][2], const T W[3][2], 
    const long long factor = 1000000000L, 
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  // T Q[3] = {0}, P[3][3] = {0}, QP[3][3] = {0};
  // characteristic_polynomials_parallel_vector2_simplex2(V, W, Q, P);
  T QP[3][3] = {0};

  std::map<long long, T> quantized_roots;
  disjoint_intervals<long long> I;
  I.set_to_complete();

  for (int i = 0; i < 3; i ++) {
    polynomial_subtraction(Q, 2, P[i], 2, QP[i]);

    const auto r0 = solve_quadratic_rational_inequality_quantized(P[i], Q, factor, epsilon);
    // std::cerr << "I: " << I << std::endl;
    // fprintf(stderr, "P: %f, %f, %f\n", P[i][0], P[i][1], P[i][2]);
    // std::cerr << "J: " << std::get<0>(r0) << std::endl;
    I.intersect(std::get<0>(r0));
    // std::cerr << "I: " << I << std::endl;
    quantized_roots.insert(std::get<1>(r0).begin(), std::get<1>(r0).end());

    const auto r1 = solve_quadratic_rational_inequality_quantized(QP[i], Q, factor, epsilon);
    // std::cerr << "I: " << I << std::endl;
    // fprintf(stderr, "P: %f, %f, %f\n", QP[i][0], QP[i][1], QP[i][2]);
    // std::cerr << "J: " << std::get<0>(r1) << std::endl;
    I.intersect(std::get<0>(r1));
    // std::cerr << "I: " << I << std::endl;
    quantized_roots.insert(std::get<1>(r1).begin(), std::get<1>(r1).end());
  }

  return std::make_tuple(I, quantized_roots);
}

template <typename T>
disjoint_intervals<T> solve_pv_inequalities_s2v2(
    const T Q[3], const T P[3][3], 
    // const T V[3][2], const T W[3][2], 
    const long long factor = 1000000000L, 
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  // const auto [I, R] = solve_pv_inequalities_quantized_s2v2(Q, P, factor, epsilon); // require c++17
  const auto tuple = solve_pv_inequalities_quantized_s2v2(Q, P, factor, epsilon);
  const auto I = std::get<0>(tuple);
  const auto R = std::get<1>(tuple); 
  return disjoint_intervals<T>(I, factor);
}

}
#endif
