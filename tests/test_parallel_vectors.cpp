#include <gtest/gtest.h>
#include <ftk/numeric/parallel_vector_solver2.hh>
#include <ftk/numeric/parallel_vector_solver3.hh>
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
