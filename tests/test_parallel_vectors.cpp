#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <ftk/numeric/parallel_vector_solver2.hh>
#include <ftk/numeric/parallel_vector_solver3.hh>
#include <ftk/numeric/rational_function.hh>
#include <ftk/numeric/rand.hh>

const int nruns = 1; // 00000;
const double epsilon = 1e-4;

TEST_CASE("solve_parallel_vector2_simplex1") {
  double V[2][2], W[2][2], lambda[2], mu[2][2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2x2(V);
    ftk::rand2x2(W);

    int nroots = ftk::solve_pv_s1v2(V, W, lambda, mu);
    for (int i = 0; i < nroots; i ++) {
      REQUIRE( ftk::verify_pv_s1v2(V, W, mu[i], epsilon) );
    }
  }
}

TEST_CASE("characteristic_polynomials_parallel_vector2_simplex2_const_w") {
  double V[3][2], w[2], P[3][2];
  double lambda, mu[3], v[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x2(V);
    ftk::rand2(w);

    ftk::characteristic_polynomials_pv_s2v2(V, w, P);

    lambda = ((double)rand() / RAND_MAX - 0.5) * 10000;
    mu[0] = ftk::polynomial_evaluate(P[0], 1, lambda); 
    mu[1] = ftk::polynomial_evaluate(P[1], 1, lambda);
    mu[2] = ftk::polynomial_evaluate(P[2], 1, lambda); // 1 - mu[0] - mu[1];

    bool succ = ftk::verify_pv_s2v2(V, w, mu, epsilon);
    if (!succ) {
      ftk::print3x2("V", V);
      fprintf(stderr, "w={%f, %f}\n", w[0], w[1]);
      fprintf(stderr, "lambda=%f, mu={%f, %f, %f}\n", lambda, mu[0], mu[1], mu[2]);
    }
    REQUIRE(succ);
  }
}

TEST_CASE("parallel_vector2_simplex2_inequality_const_w") {
  double V[3][2], w[2], P[3][2];
  double lambda, mu[3], v[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x2(V);
    ftk::rand2(w);

    ftk::characteristic_polynomials_pv_s2v2(V, w, P);
    auto I = ftk::solve_pv_inequalities_s2v2(P, w);
    if (I.empty()) continue;

    for (int k = 0; k < 3; k ++) {
      if (k == 0) lambda = I.sample();
      else if (k == 1) lambda = I.subintervals().begin()->lower();
      else if (k == 2) lambda = I.subintervals().begin()->upper();

      for (int i = 0; i < 3; i ++) 
        mu[i] = ftk::polynomial_evaluate(P[i], 1, lambda); 

      const double sum = std::abs(mu[0]) + std::abs(mu[1]) + std::abs(mu[2]);
      bool inside = (sum <= 1 + epsilon);
      if (!inside) {
        fprintf(stderr, "lambda=%f, mu=%f, %f, %f\n", lambda, mu[0], mu[1], mu[2]);
      }

      bool succ = ftk::verify_pv_s2v2(V, w, mu, epsilon);
      REQUIRE((inside && succ));
    }
  }
}

TEST_CASE("characteristic_polynomials_parallel_vector2_simplex2") {
  double V[3][2], W[3][2], Q[3], P[3][3];
  double lambda, mu[3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x2(V);
    ftk::rand3x2(W);

    ftk::characteristic_polynomials_pv_s2v2(V, W, Q, P);

    lambda = ((double)rand() / RAND_MAX - 0.5) * 10000;
    for (int i = 0; i < 3; i ++)
      mu[i] = ftk::evaluate_rational(P[i], Q, 2, lambda); 

    bool succ = ftk::verify_pv_s2v2(V, W, mu, epsilon);
    if (!succ) {
      ftk::print3x2("V", V);
      ftk::print3x2("W", W);
      ftk::print3x3("P", P);
      fprintf(stderr, "Q=%f, %f, %f\n", Q[0], Q[1], Q[2]);
      fprintf(stderr, "mu=%f, %f, %f\n", mu[0], mu[1], mu[2]);
    }
    REQUIRE(succ);
  }
}

#if 0
TEST_CASE("parallel_vector2_simplex2_inequality") {
  double V[3][2], W[3][2], Q[3], P[3][3];
  double lambda, mu[3], v[2], w[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x2(V);
    ftk::rand3x2(W);

    ftk::characteristic_polynomials_pv_s2v2(V, W, Q, P);
    auto [I, R] = ftk::solve_pv_inequalities_quantized_s2v2(Q, P);
    if (I.empty()) continue;

    std::vector<long long> quantized_lambdas;
    std::vector<double> lambdas;
    for (const auto ii : I.subintervals()) {
      lambdas.push_back(ii.lower());
      lambdas.push_back(ii.upper());
      lambdas.push_back(ii.sample());
    }

    for (int i = 0; i < lambdas.size(); i ++) {
      const double lambda = lambdas[i];
      if (!std::isnormal(lambda)) continue;

      for (int i = 0; i < 3; i ++) 
        mu[i] = ftk::polynomial_evaluate(P[i], 1, lambda); 

      const double sum = std::abs(mu[0]) + std::abs(mu[1]) + std::abs(mu[2]);
      bool inside = (sum <= 1 + epsilon);
      if (!inside) {
        fprintf(stderr, "lambda=%f, mu=%f, %f, %f\n", lambda, mu[0], mu[1], mu[2]);
      }

      bool succ = ftk::verify_pv_s2v2(V, W, mu, epsilon);
      REQUIRE(inside && succ);
    }
  }
}
#endif

TEST_CASE("characteristic_polynomials_parallel_vector3_simplex2")
{
  double V[4][3], W[4][3], Q[4], P[4][4];
  double lambda, mu[4];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand4x3(V);
    ftk::rand4x3(W);

    ftk::characteristic_polynomials_pv_s3v3(V, W, Q, P);
    lambda = ((double)rand() / RAND_MAX - 0.5) * 10000;
    mu[0] = ftk::evaluate_rational(P[0], Q, 3, lambda); 
    mu[1] = ftk::evaluate_rational(P[1], Q, 3, lambda);
    mu[2] = ftk::evaluate_rational(P[2], Q, 3, lambda); 
    mu[3] = ftk::evaluate_rational(P[3], Q, 3, lambda); 
    
    bool succ = ftk::verify_pv_s3v3(V, W, mu, epsilon);

    REQUIRE(succ);
  }
}

TEST_CASE("solve_parallel_vector3_simplex1") {
  double V[2][3], W[2][3], lambda[2], mu[2][2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2x3(V);
    ftk::rand2x3(W);

    int nroots = ftk::solve_pv_s1v3(V, W, lambda, mu);
    for (int i = 0; i < nroots; i ++) {
      REQUIRE( ftk::verify_pv_s1v3(V, W, mu[i], epsilon) );
    }
  }
}

TEST_CASE("solve_parallel_vector3_simplex2") {
  double V[3][3], W[3][3], lambda[3], mu[3][3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x3(V);
    ftk::rand3x3(W);

    int nroots = ftk::solve_pv_s2v3(V, W, lambda, mu);
    for (int i = 0; i < nroots; i ++) {
      REQUIRE( ftk::verify_pv_s2v3(V, W, mu[i], epsilon) );
    }
  }
}

int main(int argc, char **argv)
{
  Catch::Session session;
  return session.run(argc, argv);
}
