#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include <ftk/numeric/quadratic_interpolation.hh>
#include <ftk/numeric/rand.hh>
#include <ftk/external/diy/mpi.hpp>

const int nruns = 1; // 00000;
const double epsilon = 1e-9;

TEST_CASE("quadratic_interpolation") {
  double x[6][2] = {{0, 0}, {1, 0}, {1, 1}, {0.5, 0}, {1, 0.5}, {0.5, 0.5}};
  double f[6] = {0};
  double Q[3][3] = {0};
  for (int run = 0; run < nruns; run ++) {
    ftk::rand<double, 6>(f);
    ftk::quadratic_interpolation_coefficients(f, x, Q);
    for(int n=0; n<6; n++){
      double tmp_f = 0;
      for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
          double x_i = (i==2) ? 1 : x[n][i];
          double x_j = (j==2) ? 1 : x[n][j];
          tmp_f += Q[i][j] * x_i * x_j;
        }
      } 
      REQUIRE(f[n] == Catch::Approx(tmp_f).margin(epsilon));
    }
  }
}

#include "main.hh"
