#include <gtest/gtest.h>
#include <ftk/numeric/parallel_vector_solver2.hh>
#include <ftk/numeric/parallel_vector_solver3.hh>
#include <ftk/numeric/rational.hh>
#include <ftk/numeric/rand.hh>

class parallel_vectors_test : public testing::Test {
public:
  const int nruns = 100000;
  const double epsilon = 1e-4;
};

TEST_F(parallel_vectors_test, solve_parallel_vector2_simplex1) {
  double V[2][2], W[2][2], lambda[2], mu[2][2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2x2(V);
    ftk::rand2x2(W);

    int nroots = ftk::solve_parallel_vector2_simplex1(V, W, lambda, mu);
    for (int i = 0; i < nroots; i ++) {
      EXPECT_TRUE( ftk::verify_parallel_vector2_simplex1(V, W, mu[i], epsilon) );
    }
  }
}

TEST_F(parallel_vectors_test, characteristic_polynomials_parallel_vector2_simplex2_const_w) {
  double V[3][2], w[2], P[3][2];
  double lambda, mu[3], v[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x2(V);
    ftk::rand2(w);

    ftk::characteristic_polynomials_parallel_vector2_simplex2(V, w, P);

    lambda = ((double)rand() / RAND_MAX - 0.5) * 10000;
    mu[0] = ftk::polynomial_evaluate(P[0], 1, lambda); 
    mu[1] = ftk::polynomial_evaluate(P[1], 1, lambda);
    mu[2] = ftk::polynomial_evaluate(P[2], 1, lambda); // 1 - mu[0] - mu[1];

    bool succ = ftk::verify_parallel_vector2_simplex2(V, w, mu, epsilon);
    if (!succ) {
      ftk::print3x2("V", V);
      fprintf(stderr, "w={%f, %f}\n", w[0], w[1]);
      fprintf(stderr, "lambda=%f, mu={%f, %f, %f}\n", lambda, mu[0], mu[1], mu[2]);
    }
    EXPECT_TRUE(succ);
  }
}

TEST_F(parallel_vectors_test, characteristic_polynomials_parallel_vector2_simplex2) {
  double V[3][2], W[3][2], Q[3], P[3][3];
  double lambda, mu[3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x2(V);
    ftk::rand3x2(W);

    ftk::characteristic_polynomials_parallel_vector2_simplex2(V, W, Q, P);

    lambda = ((double)rand() / RAND_MAX - 0.5) * 10000;
    mu[0] = ftk::evaluate_rational(P[0], Q, 2, lambda); 
    mu[1] = ftk::evaluate_rational(P[1], Q, 2, lambda);
    mu[2] = ftk::evaluate_rational(P[2], Q, 2, lambda); // 1 - mu[0] - mu[1];

    bool succ = ftk::verify_parallel_vector2_simplex2(V, W, mu, epsilon);
    if (!succ) {
      ftk::print3x2("V", V);
      ftk::print3x2("W", W);
      ftk::print3x3("P", P);
      fprintf(stderr, "Q=%f, %f, %f\n", Q[0], Q[1], Q[2]);
      fprintf(stderr, "mu=%f, %f, %f\n", mu[0], mu[1], mu[2]);
    }
    EXPECT_TRUE(succ);
  }
}

TEST_F(parallel_vectors_test, characteristic_polynomials_parallel_vector3_simplex2)
{
  double V[4][3], W[4][3], Q[4], P[4][4];
  double lambda, mu[4];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand4x3(V);
    ftk::rand4x3(W);

    ftk::characteristic_polynomials_parallel_vector3_simplex3(V, W, Q, P);
    lambda = ((double)rand() / RAND_MAX - 0.5) * 10000;
    mu[0] = ftk::evaluate_rational(P[0], Q, 3, lambda); 
    mu[1] = ftk::evaluate_rational(P[1], Q, 3, lambda);
    mu[2] = ftk::evaluate_rational(P[2], Q, 3, lambda); 
    mu[3] = ftk::evaluate_rational(P[3], Q, 3, lambda); 
    
    bool succ = ftk::verify_parallel_vector3_simplex3(V, W, mu, epsilon);

    EXPECT_TRUE(succ);
  }
}

TEST_F(parallel_vectors_test, solve_parallel_vector3_simplex1) {
  double V[2][3], W[2][3], lambda[2], mu[2][2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2x3(V);
    ftk::rand2x3(W);

    int nroots = ftk::solve_parallel_vector3_simplex1(V, W, lambda, mu);
    for (int i = 0; i < nroots; i ++) {
      EXPECT_TRUE( ftk::verify_parallel_vector3_simplex1(V, W, mu[i], epsilon) );
    }
  }
}

TEST_F(parallel_vectors_test, solve_parallel_vector3_simplex2) {
  double V[3][3], W[3][3], lambda[3], mu[3][3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x3(V);
    ftk::rand3x3(W);

    int nroots = ftk::solve_parallel_vector3_simplex2(V, W, lambda, mu);
    for (int i = 0; i < nroots; i ++) {
      EXPECT_TRUE( ftk::verify_parallel_vector3_simplex2(V, W, mu[i], epsilon) );
    }
  }
}
