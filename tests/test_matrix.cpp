#include <gtest/gtest.h>
#include <ftk/numeric/matrix_inverse.hh>
#include <ftk/numeric/matrix_multiplication.hh>
#include <ftk/numeric/eigen_solver2.hh>
#include <ftk/numeric/eigen_solver3.hh>
#include <ftk/numeric/rand.hh>
#include <ftk/numeric/linear_solver.hh>
#include <ftk/numeric/print.hh>

class matrix_test : public testing::Test {
public:
  const int nruns = 100000;
  const double epsilon = 1e-4;
};

TEST_F(matrix_test, matrix_inverse2) {
  double A[2][2], invA[2][2], I[2][2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2x2(A);
    ftk::matrix_inverse2x2(A, invA);
    ftk::matrix2x2_matrix2x2_multiplication(A, invA, I);

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
    ftk::matrix3x3_matrix3x3_multiplication(A, invA, I);

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

TEST_F(matrix_test, solve_linear2x2) {
  double A[2][2], b[2], x[2], bb[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2x2(A);
    ftk::rand2(b);
    ftk::solve_linear2x2(A, b, x);

    ftk::matrix2x2_vector2_multiplication(A, x, bb);

    EXPECT_NEAR(b[0], bb[0], epsilon);
    EXPECT_NEAR(b[1], bb[1], epsilon);
  }
}

TEST_F(matrix_test, solve_linear3x3) {
  double A[3][3], b[3], x[3], bb[3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x3(A);
    ftk::rand3(b);
    ftk::solve_linear3x3(A, b, x);

    ftk::matrix3x3_vector3_multiplication(A, x, bb);

    EXPECT_NEAR(b[0], bb[0], epsilon);
    EXPECT_NEAR(b[1], bb[1], epsilon);
    EXPECT_NEAR(b[2], bb[2], epsilon);
  }
}

TEST_F(matrix_test, solve_eigenvalues_symmetric2x2) {
  double A[2][2], eig[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand_symmetric2x2(A);
    ftk::solve_eigenvalues_symmetric2x2(A, eig);

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

TEST_F(matrix_test, solve_eigenvalues2x2) {
  double A[2][2], eig[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2x2(A);
    int nroots = ftk::solve_eigenvalues2x2(A, eig);

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

TEST_F(matrix_test, solve_generalized_eigenvalues2x2) {
  double A[2][2], B[2][2], eig[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2x2(A); 
    ftk::rand2x2(B);
    const int n = ftk::solve_generalized_eigenvalues2x2(A, B, eig);

    for (int i = 0; i < n; i ++) {
      double M[2][2];
      for (int j = 0; j < 2; j ++) 
        for (int k = 0; k < 2; k ++) 
          M[j][k] = A[j][k] - eig[i] * B[j][k];

      const double det = ftk::det2(M);
      EXPECT_NEAR(0.0, det, epsilon);
    }
  }
}

TEST_F(matrix_test, solve_eigenvectors2x2)
{
  double A[2][2], eig[2], eigvecs[2][2], x[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2x2(A);
    int nroots = ftk::solve_eigenvalues2x2(A, eig);
    ftk::solve_eigenvectors2x2(A, nroots, eig, eigvecs);

    // ftk::print2x2("A", A);
    for (int i = 0; i < nroots; i ++) {
      ftk::matrix2x2_vector2_multiplication(A, eigvecs[i], x);
      // fprintf(stderr, "i=%d, eig=%f, eigvec={%f, %f}, x={%f, %f}\n", 
      //     i, eig[i], eigvecs[i][0], eigvecs[i][1], x[0], x[1]);

      EXPECT_NEAR(eig[i] * eigvecs[i][0], x[0], epsilon);
      EXPECT_NEAR(eig[i] * eigvecs[i][1], x[1], epsilon);
    }
  }
}

TEST_F(matrix_test, solve_eigenvalues_symmetric3x3) {
  double A[3][3], eig[3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand_symmetric3x3(A);
    ftk::solve_eigenvalues_symmetric3x3(A, eig);

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

TEST_F(matrix_test, solve_eigenvectors3x3)
{
  const int nroots = 3; // for now, we use symetric matrix for testing
  double A[3][3], eig[3], eigvecs[3][3], x[3], c[3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand_symmetric3x3(A);
    ftk::solve_eigenvalues_symmetric3x3(A, eig);
    ftk::solve_eigenvectors3x3(A, nroots, eig, eigvecs);

    for (int i = 0; i < nroots; i ++) {
      ftk::matrix3x3_vector3_multiplication(A, eigvecs[i], x);

      EXPECT_NEAR(eig[i] * eigvecs[i][0], x[0], epsilon);
      EXPECT_NEAR(eig[i] * eigvecs[i][1], x[1], epsilon);
      EXPECT_NEAR(eig[i] * eigvecs[i][2], x[2], epsilon);
    }
  }
}

TEST_F(matrix_test, solve_generalized_eigenvalues3x3) {
  double A[3][3], B[3][3], eig[3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x3(A);
    ftk::rand3x3(B);
    const int n = ftk::solve_generalized_eigenvalues3x3(A, B, eig);

    for (int i = 0; i < n; i ++) {
      double M[3][3];
      for (int j = 0; j < 3; j ++) 
        for (int k = 0; k < 3; k ++) 
          M[j][k] = A[j][k] - eig[i] * B[j][k];

      const double det = ftk::det3(M);
      EXPECT_NEAR(0.0, det, epsilon);
    }
  }
}

