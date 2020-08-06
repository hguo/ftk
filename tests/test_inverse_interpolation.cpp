#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/rand.hh>
#include <ftk/numeric/trilinear_interpolation.hh>
#include <ftk/numeric/inverse_trilinear_interpolation_solver.hh>

const int nruns = 1; // 00000;
const double epsilon = 1e-9;

TEST_CASE("inverse_linear_interpolation_2simplex_vector2") {
  double V[3][2], mu[3], v[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x2(V);
    ftk::inverse_lerp_s2v2(V, mu);
    ftk::lerp_s2v2(V, mu, v);

    REQUIRE(v[0] == Approx(0.0).margin(epsilon));
    REQUIRE(v[1] == Approx(0.0).margin(epsilon));
  }
}

TEST_CASE("inverse_linear_interpolation_3simplex_vector3") {
  double V[4][3], mu[4], v[3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand4x3(V);
    ftk::inverse_lerp_s3v3(V, mu);
    ftk::lerp_s3v3(V, mu, v);

    REQUIRE(v[0] == Approx(0.0).margin(epsilon));
    REQUIRE(v[1] == Approx(0.0).margin(epsilon));
    REQUIRE(v[2] == Approx(0.0).margin(epsilon));
  }
}

TEST_CASE("trilinear_interpolation3") {
  double V[8][3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand<double, 8, 3>(V);
    for(int i=0; i<8; i++){
      double alpha = (i>>2) & 1;
      double beta = (i>>1) & 1;
      double gamma = i & 1;
      double result[3];
      ftk::trilinear_interpolation3(V, alpha, beta, gamma, result);
      REQUIRE(V[i][0] == Approx(result[0]).margin(epsilon));
      REQUIRE(V[i][1] == Approx(result[1]).margin(epsilon));
      REQUIRE(V[i][2] == Approx(result[2]).margin(epsilon));
    }
  }
}

#if FTK_HAVE_MPSOLVE
TEST_CASE("inverse_trilinear_interpolation3") {
  double V[8][3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand<double, 8, 3>(V);
    double result[6][3];
    int num = ftk::inverse_trilinear_interpolation8_3(V, result);
    for(int i=0; i<num; i++){
      double value[3];
      ftk::trilinear_interpolation3(V, result[i][0], result[i][1], result[i][2], value);
      REQUIRE(value[0] == Approx(0.0).margin(epsilon));
      REQUIRE(value[1] == Approx(0.0).margin(epsilon));
      REQUIRE(value[2] == Approx(0.0).margin(epsilon));
    }
  }
}
#endif


int main(int argc, char **argv)
{
  Catch::Session session;
  return session.run(argc, argv);
}
