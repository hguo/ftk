#include <gtest/gtest.h>
#include <ftk/numeric/matrix_inverse.hh>
#include <ftk/numeric/matrix_multiplication.hh>
#include <ftk/numeric/eigen_solver.hh>
#include <ftk/numeric/rand.hh>
#include <ftk/numeric/print.hh>

class matrix_test : public testing::Test {
public:
  const int nruns = 100000;
  const double epsilon = 1e-10;
};

TEST_F(matrix_test, matrix_inverse2) {
  double A[2][2], invA[2][2], I[2][2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2x2(A);
    ftk::matrix_inverse2x2(A, invA);
    ftk::matrix_matrix_multiplication_2x2_2x2(A, invA, I);

    EXPECT_NEAR(1.0, I[0][0], epsilon);
    EXPECT_NEAR(0.0, I[0][1], epsilon);
    EXPECT_NEAR(0.0, I[1][0], epsilon);
    EXPECT_NEAR(1.0, I[1][1], epsilon);
  }
}

TEST_F(matrix_test, matrix_inverse3) {
  double A[3][3], invA[3][3], I[3][3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x3(A);
    ftk::matrix_inverse3x3(A, invA);
    ftk::matrix_matrix_multiplication_3x3_3x3(A, invA, I);

    EXPECT_NEAR(1.0, I[0][0], epsilon);
    EXPECT_NEAR(0.0, I[0][1], epsilon);
    EXPECT_NEAR(0.0, I[0][2], epsilon);
    EXPECT_NEAR(0.0, I[1][0], epsilon);
    EXPECT_NEAR(1.0, I[1][1], epsilon);
    EXPECT_NEAR(0.0, I[1][2], epsilon);
    EXPECT_NEAR(0.0, I[2][0], epsilon);
    EXPECT_NEAR(0.0, I[2][1], epsilon);
    EXPECT_NEAR(1.0, I[2][2], epsilon);
  }
}

TEST_F(matrix_test, solve_eigenvalues_real_symmetric2x2) {
  double A[2][2], eig[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand_symmetric2x2(A);
    ftk::solve_eigenvalues_real_symmetric2x2(A, eig);

    for (int i = 0; i < 2; i ++) {
      double B[2][2] = {
        {A[0][0] - eig[i], A[0][1]}, 
        {A[1][0], A[1][1] - eig[i]}
      };

      const double det = ftk::det2(B);
      EXPECT_NEAR(0.0, det, epsilon);
    }
  }
}

TEST_F(matrix_test, solve_real_eigenvalues_real2x2) {
  double A[2][2], eig[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2x2(A);
    int nroots = ftk::solve_real_eigenvalues_real2x2(A, eig);

    for (int i = 0; i < nroots; i ++) {
      double B[2][2] = {
        {A[0][0] - eig[i], A[0][1]}, 
        {A[1][0], A[1][1] - eig[i]}
      };

      const double det = ftk::det2(B);
      EXPECT_NEAR(0.0, det, epsilon);
    }
  }
}

TEST_F(matrix_test, solve_eigenvalues_real_symmetric3x3) {
  double A[3][3], eig[3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand_symmetric3x3(A);
    ftk::solve_eigenvalues_real_symmetric3x3(A, eig);

    for (int i = 0; i < 3; i ++) {
      double B[3][3] = {
        {A[0][0] - eig[i], A[0][1], A[0][2]}, 
        {A[1][0], A[1][1] - eig[i], A[1][2]},
        {A[2][0], A[2][1], A[2][2] - eig[i]}
      };

      const double det = ftk::det3(B);
      EXPECT_NEAR(0.0, det, epsilon);
    }
  }
}
