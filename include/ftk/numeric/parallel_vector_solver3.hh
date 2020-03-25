#ifndef _FTK_PARALLEL_VECTOR_SOLVER3_H
#define _FTK_PARALLEL_VECTOR_SOLVER3_H

#include <ftk/numeric/matrix_inverse.hh>
#include <ftk/numeric/matrix_addition.hh>
#include <ftk/numeric/eigen_solver2.hh>
#include <ftk/numeric/eigen_solver3.hh>
#include <ftk/numeric/linear_solver.hh>
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
inline bool verify_pv_s0v3(const T v[3], const T w[3], const T epsilon = std::numeric_limits<T>::epsilon())
{
  double c[3];
  ftk::cross_product(v, w, c);
  
  const double norm = ftk::vector_2norm_3(c);
  return norm <= epsilon;
}

template <typename T>
inline bool verify_pv_s1v3(
    const T V[2][3], const T W[2][3], const T mu[2], 
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  double v[3], w[3];
  ftk::lerp_s1v3(V, mu, v);
  ftk::lerp_s1v3(W, mu, w);

  return verify_pv_s0v3(v, w, epsilon);
}

template <typename T>
inline bool verify_pv_s2v3(const T V[3][3], const T W[3][3], const T mu[3], 
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  double v[3], w[3];
  ftk::lerp_s2v3(V, mu, v);
  ftk::lerp_s2v3(W, mu, w);

  return verify_pv_s0v3(v, w, epsilon);
}

template <typename T>
inline bool verify_pv_s3v3(const T V[4][3], const T W[4][3], const T mu[4],
    const T epsilon = std::numeric_limits<T>::epsilon())
{
  double v[3], w[3];
  ftk::lerp_s3v3(V, mu, v);
  ftk::lerp_s3v3(W, mu, w);

  return verify_pv_s0v3(v, w, epsilon);
}

template <typename T>
inline int solve_pv_s2v3(const T VV[3][3], const T WW[3][3], 
    T lambda[3], T mu[3][3], const T epsilon = std::numeric_limits<T>::epsilon()) // was 1e-6
{
  T V[3][3], W[3][3]; // transposed V and W
  transpose3x3(VV, V);
  transpose3x3(WW, W);

  T P[4]; // characteristic polynomial
  characteristic_polynomial_3x3(V, W, P);

  T l[3] = {0};
  const int n_roots = solve_cubic_real(P, l, epsilon);
  int n_solutions = 0;

  for (int i = 0; i < n_roots; i ++) {
    if (std::abs(l[i]) <= epsilon) continue;

    const T M[3][2] = {
      {(V[0][0] - V[0][2]) - l[i] * (W[0][0] - W[0][2]), (V[0][1] - V[0][2]) - l[i] * (W[0][1] - W[0][2])}, 
      {(V[1][0] - V[1][2]) - l[i] * (W[1][0] - W[1][2]), (V[1][1] - V[1][2]) - l[i] * (W[1][1] - W[1][2])},
      {(V[2][0] - V[2][2]) - l[i] * (W[2][0] - W[2][2]), (V[2][1] - V[2][2]) - l[i] * (W[2][1] - W[2][2])}
    };
    const T b[3] = {
      -(V[0][2] - l[i]*W[0][2]), 
      -(V[1][2] - l[i]*W[1][2]),
      -(V[2][2] - l[i]*W[2][2])
    };
    T nu[3];
    // const auto det = solve_linear_real2x2(M, b, nu, T(0)); // epsilon);
    const auto cond = solve_least_square3x2(M, b, nu, T(0));
    nu[2] = T(1) - nu[0] - nu[1];
    
    // fprintf(stderr, "lambda=%f, nu={%.20f, %.20f, %.20f}, cond=%.20f\n", l[i], nu[0], nu[1], nu[2], cond);
    // print3x2("M", M);
    // fprintf(stderr, "b={%f, %f, %f}\n", b[0], b[1], b[2]);

    // if (det == T(0)) {
    //   fprintf(stderr, "FATAL: det(M)=0, P={%f, %f, %f, %f}, n_roots=%d, i=%d, lambda=%f, nu={%f, %f, %f}\n", 
    //       P[0], P[1], P[2], P[3], n_roots, i, l[i], nu[0], nu[1], nu[2]);
    //   print3x3("V", VV);
    //   print3x3("W", WW);
    // }
    // if (std::abs(det) <= epsilon) continue; // FIXME: fixed threshold
    if (nu[0] >= -epsilon && nu[0] <= 1+epsilon && nu[1] >= -epsilon && nu[1] <= 1+epsilon && nu[2] >= -epsilon && nu[2] <= 1+epsilon) {
    // if (nu[0] >= 0 && nu[0] <= 1 && nu[1] >= 0 && nu[1] <= 1 && nu[2] >= 0 && nu[2] <= 1) {
#if 0
      if (cond > 1e9) { // bad condition number
        double v[3], w[3], c[3]; 
        ftk::linear_interpolation3_3(VV, nu, v);
        ftk::linear_interpolation3_3(WW, nu, w);
        ftk::cross_product(v, w, c);
        const double norm = ftk::vector_2norm_3(c);
#if 0
        fprintf(stderr, "rejecting by condition number, P={%f, %f, %f, %f}, n_roots=%d, roots={%f, %f, %f}, lambda=%f, nu={%f, %f, %f}, cond=%f, ||v x w||=%f\n", 
            P[0], P[1], P[2], P[3], n_roots, l[0], l[1], l[2], 
            l[i], nu[0], nu[1], nu[2], cond, norm);
        print3x3("V", VV);
        print3x3("W", WW);
        print3x2("M", M);
        print3("b", b);
#endif
        continue; // bad condition number
      }
#endif

#if 1
      // check interpolation results
      double v[3], w[3], c[3]; 
      ftk::lerp_s2v3(VV, nu, v);
      ftk::lerp_s2v3(WW, nu, w);
      ftk::cross_product(v, w, c);
      const double norm = ftk::vector_2norm_3(c);
      if (norm > 1e-2) {
        fprintf(stderr, "rejecting: nu={%f, %f, %f}, v={%f, %f, %f}, w={%f, %f, %f}, c={%f, %f, %f}, norm=%f\n", 
            nu[0], nu[1], nu[2], v[0], v[1], v[2], w[0], w[1], w[2], c[0], c[1], c[2], norm);
        print3x3("V", VV);
        print3x3("W", WW);
        continue; // reject
      }
#endif
      
      const int j = n_solutions ++;
      lambda[j] = l[i];
      mu[j][0] = nu[0]; 
      mu[j][1] = nu[1];
      mu[j][2] = nu[2];
    }
  }
 
  // if (n_solutions > 0) 
  //   if (P[3] == T(0) || std::isinf(P[3]) || std::isnan(P[3])) {
  //     fprintf(stderr, "WARN: degeneracy case, P={%f, %f, %f, %f}, n_roots=%d, n_solutions=%d\n", 
  //         P[0], P[1], P[2], P[3], n_roots, n_solutions); 
  //   }

  return n_solutions;
}

template <typename T>
inline int solve_pv_s1v3(
    const T V[2][3], const T W[2][3], 
    T lambda[2], T mu[2][2], 
    const T epsilon = std::numeric_limits<T>::epsilon()) // 1e-9)
{
  const T P[3][2] = { // rhs
    {V[1][0], -W[1][0]},
    {V[1][1], -W[1][1]},
    {V[1][2], -W[1][2]}
  };
  const T norm = matrix_frobenius_norm_real3x2(P);
  if (norm <= epsilon)
    return std::numeric_limits<int>::max(); // unlimited.

  const T Q[3][2] = { // lhs
    {V[1][0]-V[0][0], W[0][0]-W[1][0]},
    {V[1][1]-V[0][1], W[0][1]-W[1][1]},
    {V[1][2]-V[0][2], W[0][2]-W[1][2]}
  };

  // print2x3("V", V);
  // print2x3("W", W);
  //print3x2("P", P);
  //print3x2("Q", Q);

  T E[3][3], F[3][3];
  polynomial_multiplication(P[1], 1, Q[2], 1, E[0]);
  polynomial_multiplication(P[2], 1, Q[1], 1, F[0]);
  polynomial_multiplication(P[0], 1, Q[2], 1, E[1]);
  polynomial_multiplication(P[2], 1, Q[0], 1, F[1]);
  polynomial_multiplication(P[0], 1, Q[1], 1, E[2]);
  polynomial_multiplication(P[1], 1, Q[0], 1, F[2]);

  T H[3][3];
  for (int i = 0; i < 3; i ++)
    polynomial_subtraction(E[i], 2, F[i], 2, H[i]);

  T C[3];
  cross_product(H[0], H[1], C);
  // fprintf(stderr, "E={%f, %f, %f}, F={%f, %f, %f}, C={%f, %f, %f}, normC=%.20f\n", 
  //     E[0], E[1], E[2], F[0], F[1], F[2], C[0], C[1], C[2], vector_2norm_3(C));

  if (vector_2norm_3(H[0]) > epsilon && vector_2norm_3(H[1]) > epsilon && vector_2norm_3(C) <= epsilon) {
    T l[2], nu[2];
    const int n_roots = solve_quadratic_real(H[0], l);
    // fprintf(stderr, "n_roots=%d, l={%f, %f}\n", n_roots, l[0], l[1]);
    int n = 0;
    for (int i = 0; i < n_roots; i ++) {
      if (std::abs(l[i]) <= epsilon) continue;
      // if (std::abs(polynomial_evaluate(Q[0], 1, l[i])) < epsilon) continue;
#if 0
      const T A[3] = {
        polynomial_evaluate(P[0], 1, l[i]), 
        polynomial_evaluate(P[1], 1, l[i]), 
        polynomial_evaluate(P[2], 1, l[i])
      }, B[3] = {
        polynomial_evaluate(Q[0], 1, l[i]), 
        polynomial_evaluate(Q[1], 1, l[i]), 
        polynomial_evaluate(Q[2], 1, l[i])
      };
      fprintf(stderr, "A={%f, %f, %f}, B={%f, %f, %f}\n", 
          A[0], A[1], A[2], B[0], B[1], B[2]);
#endif
      nu[i] = polynomial_evaluate(P[0], 1, l[i]) / polynomial_evaluate(Q[0], 1, l[i]);
      // fprintf(stderr, "l=%f, nu=%.20f\n", l[i], nu[i]);
      if (nu[i] >= -epsilon && nu[i] <= 1+epsilon) {
        lambda[n] = l[i];
        mu[n][0] = nu[i];
        mu[n][1] = 1 - nu[i];
        n ++;
      }
    }
    return n;
  } else 
    return 0;
}

template <typename T>
inline void characteristic_polynomials_pv_s3v3(const T V[4][3], const T W[4][3], T Q[4], T P[4][4])
{
  const T A[3][3] = { // linear transformation
    {V[0][0] - V[3][0], V[1][0] - V[3][0], V[2][0] - V[3][0]}, 
    {V[0][1] - V[3][1], V[1][1] - V[3][1], V[2][1] - V[3][1]}, 
    {V[0][2] - V[3][2], V[1][2] - V[3][2], V[2][2] - V[3][2]}
  };
  const T B[3][3] = {
    {W[0][0] - W[3][0], W[1][0] - W[3][0], W[2][0] - W[3][0]}, 
    {W[0][1] - W[3][1], W[1][1] - W[3][1], W[2][1] - W[3][1]}, 
    {W[0][2] - W[3][2], W[1][2] - W[3][2], W[2][2] - W[3][2]}
  };
  const T a[3] = {V[3][0], V[3][1], V[3][2]};
  const T b[3] = {W[3][0], W[3][1], W[3][2]};

  const T rhs[3][2] = { // polynomials for the right hand side vectors
    {-a[0], b[0]}, 
    {-a[1], b[1]}, 
    {-a[2], b[2]}};

  // coefficients for Q
  characteristic_polynomial_3x3(A, B, Q);

  T adj[3][3][3]; // build adjugate matrix
  characteristic_polynomial_2x2(A[1][1], A[1][2], A[2][1], A[2][2], B[1][1], B[1][2], B[2][1], B[2][2], adj[0][0]);
  characteristic_polynomial_2x2(A[1][0], A[1][2], A[2][0], A[2][2], B[1][0], B[1][2], B[2][0], B[2][2], adj[1][0]);
  characteristic_polynomial_2x2(A[1][0], A[1][1], A[2][0], A[2][1], B[1][0], B[1][1], B[2][0], B[2][1], adj[2][0]);
  characteristic_polynomial_2x2(A[0][1], A[0][2], A[2][1], A[2][2], B[0][1], B[0][2], B[2][1], B[2][2], adj[0][1]);
  characteristic_polynomial_2x2(A[0][0], A[0][2], A[2][0], A[2][2], B[0][0], B[0][2], B[2][0], B[2][2], adj[1][1]);
  characteristic_polynomial_2x2(A[0][0], A[0][1], A[2][0], A[2][1], B[0][0], B[0][1], B[2][0], B[2][1], adj[2][1]);
  characteristic_polynomial_2x2(A[0][1], A[0][2], A[1][1], A[1][2], B[0][1], B[0][2], B[1][1], B[1][2], adj[0][2]);
  characteristic_polynomial_2x2(A[0][0], A[0][2], A[1][0], A[1][2], B[0][0], B[0][2], B[1][0], B[1][2], adj[1][2]);
  characteristic_polynomial_2x2(A[0][0], A[0][1], A[1][0], A[1][1], B[0][0], B[0][1], B[1][0], B[1][1], adj[2][2]);
  polynomial_scalar_multiplication(adj[0][1], 2, T(-1));
  polynomial_scalar_multiplication(adj[1][0], 2, T(-1));
  polynomial_scalar_multiplication(adj[1][2], 2, T(-1));
  polynomial_scalar_multiplication(adj[2][1], 2, T(-1));

  for (int i = 0; i < 3; i ++) {
    T poly[4] = {0};
    for (int j = 0; j < 4; j ++) 
      P[i][j] = 0; // clear results
    for (int j = 0; j <3; j ++) {
      polynomial_multiplication(adj[i][j], 2, rhs[j], 1, poly);
      polynomial_addition_in_place(P[i], 3, poly, 3);
    }
  }

  P[3][0] = Q[0] - P[0][0] - P[1][0] - P[2][0];
  P[3][1] = Q[1] - P[0][1] - P[1][1] - P[2][1];
  P[3][2] = Q[2] - P[0][2] - P[1][2] - P[2][2];
  P[3][3] = Q[3] - P[0][3] - P[1][3] - P[2][3];
}

template <typename T>
std::tuple<disjoint_intervals<long long>, std::map<long long, T> >
solve_pv_inequalities_quantized_s3v3(
    const T V[4][3], const T W[4][3], const long long factor = 1000000000L, 
    const T epsilon=std::numeric_limits<T>::epsilon())
{
  T Q[4] = {0}, P[4][4] = {0}, QP[4][4] = {0};
  characteristic_polynomials_pv_s3v3(V, W, Q, P);

  std::map<long long, T> quantized_roots;
  disjoint_intervals<long long> I;
  I.set_to_complete();

  for (int i = 0; i < 4; i ++) {
    polynomial_subtraction(Q, 3, P[i], 3, QP[i]); // QP[i] = Q - P[i]

    // fprintf(stderr, "Q={%.20f, %.20f, %.20f, %.20f}, P[%d]={%.20f, %.20f, %.20f, %.20f}, QP[%d]={%.20f, %.20f, %.20f, %.20f}\n", 
    //     Q[0], Q[1], Q[2], Q[3], 
    //     i, P[i][0], P[i][1], P[i][2], P[i][3], 
    //     i, QP[i][0], QP[i][1], QP[i][2], QP[i][3]);

    const auto r0 = solve_cubic_rational_inequality_quantized(P[i], Q, factor, epsilon);
    // std::cerr << "I0: " << disjoint_intervals<T>(std::get<0>(r0), factor) << std::endl;
    // std::cerr << "I0: " << std::get<0>(r0) << std::endl;
    
    I.intersect(std::get<0>(r0));
    quantized_roots.insert(std::get<1>(r0).begin(), std::get<1>(r0).end());
    // std::cerr << "I: " << disjoint_intervals<T>(I, factor) << std::endl;
    
    const auto r1 = solve_cubic_rational_inequality_quantized(QP[i], Q, factor, epsilon);
    // std::cerr << "I1: " << disjoint_intervals<T>(std::get<0>(r1), factor) << std::endl;
    
    I.intersect(std::get<0>(r1));
    quantized_roots.insert(std::get<1>(r1).begin(), std::get<1>(r1).end());
    // std::cerr << "I: " << disjoint_intervals<T>(I, factor) << std::endl;
  }

  return std::make_tuple(I, quantized_roots);
}

template <typename T>
disjoint_intervals<T> solve_pv_inequalities_s3v3(
    const T V[4][3], const T W[4][3], const long long factor = 1000000000L, const T epsilon=std::numeric_limits<T>::epsilon())
{
  // const auto [I, R] = solve_pv_inequalities_quantized_s3v3(V, W, factor, epsilon); // require c++17
  const auto tuple = solve_pv_inequalities_quantized_s3v3(V, W, factor, epsilon);
  const auto I = std::get<0>(tuple);
  const auto R = std::get<1>(tuple);
  return disjoint_intervals<T>(I, factor);
}

}
#endif
