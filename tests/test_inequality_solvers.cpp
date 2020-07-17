#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <ftk/numeric/rand.hh>
#include <ftk/numeric/polynomial.hh>
#include <ftk/numeric/linear_inequality_solver.hh>

const int nruns = 1; // 00000;
const double epsilon = 1e-9;

TEST_CASE("linear_inequality_test") {
  double P[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2(P);
    auto I = ftk::solve_linear_inequality(P[1], P[0]);
    auto x = I.sample();
    auto y = ftk::polynomial_evaluate(P, 1, x);

    REQUIRE(y >= 0.0);
  }
}

TEST_CASE("linear_inequality_test2") {
  const int n = 2;
  double P[n][2];
  ftk::disjoint_intervals<double> I;

  for (int run = 0; run < nruns; run ++) {
    I.set_to_complete();
    ftk::rand<double, n, 2>(P);
    for (int i = 0; i < n; i ++) {
      auto ii = ftk::solve_linear_inequality(P[i][1], P[i][0]);
      I.intersect(ii);
    }

    if (I.empty()) continue; // TODO: check false positives
    // std::cerr << I << std::endl;
    
    auto x = I.sample();
    for (int i = 0; i < n; i ++) {
      auto y = ftk::polynomial_evaluate(P[i], 1, x);
      REQUIRE(y >= 0.0);
    }
  }
}

int main(int argc, char **argv)
{
  Catch::Session session;
  return session.run(argc, argv);
}
