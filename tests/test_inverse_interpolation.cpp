#include <gtest/gtest.h>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/rand.hh>
#include <ftk/numeric/trilinear_interpolation.hh>
#include <ftk/numeric/inverse_trilinear_interpolation_solver.hh>

class inverse_interpolation_test : public testing::Test {
public:
  const int nruns = 1; // 00000;
  const double epsilon = 1e-9;
};

TEST_F(inverse_interpolation_test, inverse_linear_interpolation_2simplex_vector2) {
  double V[3][2], mu[3], v[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x2(V);
    ftk::inverse_lerp_s2v2(V, mu);
    ftk::lerp_s2v2(V, mu, v);

    EXPECT_NEAR(0.0, v[0], epsilon);
    EXPECT_NEAR(0.0, v[1], epsilon);
  }
}

TEST_F(inverse_interpolation_test, inverse_linear_interpolation_3simplex_vector3) {
  double V[4][3], mu[4], v[3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand4x3(V);
    ftk::inverse_lerp_s3v3(V, mu);
    ftk::lerp_s3v3(V, mu, v);

    EXPECT_NEAR(0.0, v[0], epsilon);
    EXPECT_NEAR(0.0, v[1], epsilon);
    EXPECT_NEAR(0.0, v[2], epsilon);
  }
}

TEST_F(inverse_interpolation_test, trilinear_interpolation3) {
  double V[8][3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand<double, 8, 3>(V);
    for(int i=0; i<8; i++){
      double alpha = (i>>2) & 1;
      double beta = (i>>1) & 1;
      double gamma = i & 1;
      double result[3];
      ftk::trilinear_interpolation3(V, alpha, beta, gamma, result);
      EXPECT_NEAR(0.0, V[i][0] - result[0], epsilon);
      EXPECT_NEAR(0.0, V[i][1] - result[1], epsilon);
      EXPECT_NEAR(0.0, V[i][2] - result[2], epsilon);
    }
  }
}

TEST_F(inverse_interpolation_test, inverse_trilinear_interpolation3) {
  double V[8][3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand<double, 8, 3>(V);
    double result[6][3];
    int num = ftk::inverse_trilinear_interpolation8_3(V, result);
    for(int i=0; i<num; i++){
      double value[3];
      ftk::trilinear_interpolation3(V, result[i][0], result[i][1], result[i][2], value);
      EXPECT_NEAR(0.0, value[0], epsilon);
      EXPECT_NEAR(0.0, value[1], epsilon);
      EXPECT_NEAR(0.0, value[2], epsilon);
    }
  }
}
