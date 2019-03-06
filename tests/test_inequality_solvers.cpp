#include <gtest/gtest.h>
#include <ftk/numeric/rand.hh>
#include <ftk/numeric/polynomial.hh>
#include <ftk/numeric/linear_inequality_solver.hh>

class inequality_test : public testing::Test {
public:
  const int nruns = 100000;
  const double epsilon = 1e-9;
};

TEST_F(inequality_test, linear_inequality_test) {
  double P[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2(P);
    auto I = ftk::solve_linear_inequality(P[1], P[0]);
    auto x = I.sample();
    auto y = ftk::polynomial_evaluate(P, 1, x);

    EXPECT_GE(y, 0.0);
  }
}
