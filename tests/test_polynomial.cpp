#include <gtest/gtest.h>
#include <ftk/numeric/polynomial.hh>
#include <ftk/numeric/quadratic_solver.hh>
#include <ftk/numeric/cubic_solver.hh>
#include <random>

class polynomial_test : public testing::Test {
public:
  const int nruns = 100000;
  const double epsilon = 1e-10;
  
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> d{0, 1};
};

TEST_F(polynomial_test, solve_quadratic_real)
{
  double P[3], x[2];
  for (int run = 0; run < nruns; run ++) {
    for (int i = 0; i < 3; i ++)
      P[i] = d(gen);
    
    int nroots = ftk::solve_quadratic_real(P, x);
    for (int i = 0; i < nroots; i ++) {
      double val = ftk::polynomial_evaluate(P, 2, x[i]);
      EXPECT_NEAR(0.0, val, epsilon);
    }
  }
}

#if 0
TEST_F(polynomial_test, solve_cubic_real)
{
  double P[4], x[3];
  for (int run = 0; run < nruns; run ++) {
    for (int i = 0; i < 4; i ++)
      P[i] = d(gen);
    
    int nroots = ftk::solve_cubic_real(P, x);
    for (int i = 0; i < nroots; i ++) {
      double val = ftk::polynomial_evaluate(P, 3, x[i]);
      EXPECT_NEAR(0.0, val, epsilon);
    }
  }
}
#endif
