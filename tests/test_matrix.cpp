#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <ftk/numeric/matrix_inverse.hh>
#include <ftk/numeric/matrix_multiplication.hh>
#include <ftk/numeric/eigen_solver2.hh>
#include <ftk/numeric/eigen_solver3.hh>
#include <ftk/numeric/rand.hh>
#include <ftk/numeric/linear_solver.hh>
#include <ftk/numeric/print.hh>

const int nruns = 1000;
const double epsilon = 1e-3;

TEST_CASE("matrix_inverse2") {
  double A[2][2], invA[2][2], I[2][2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2x2(A);
    ftk::matrix_inverse2x2(A, invA);
    ftk::matrix2x2_matrix2x2_multiplication(A, invA, I);

    REQUIRE(I[0][0] == Approx(1.0).margin(epsilon));
    REQUIRE(I[0][1] == Approx(0.0).margin(epsilon));
    REQUIRE(I[1][0] == Approx(0.0).margin(epsilon));
    REQUIRE(I[1][1] == Approx(1.0).margin(epsilon));
  }
}

TEST_CASE("matrix_inverse3") {
  double A[3][3], invA[3][3], I[3][3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x3(A);
    ftk::matrix_inverse3x3(A, invA);
    ftk::matrix3x3_matrix3x3_multiplication(A, invA, I);

    REQUIRE(I[0][0] == Approx(1.0).margin(epsilon));
    REQUIRE(I[0][1] == Approx(0.0).margin(epsilon));
    REQUIRE(I[0][2] == Approx(0.0).margin(epsilon));
    REQUIRE(I[1][0] == Approx(0.0).margin(epsilon));
    REQUIRE(I[1][1] == Approx(1.0).margin(epsilon));
    REQUIRE(I[1][2] == Approx(0.0).margin(epsilon));
    REQUIRE(I[2][0] == Approx(0.0).margin(epsilon));
    REQUIRE(I[2][1] == Approx(0.0).margin(epsilon));
    REQUIRE(I[2][2] == Approx(1.0).margin(epsilon));
  }
}

TEST_CASE("solve_linear2x2") {
  double A[2][2], b[2], x[2], bb[2];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand2x2(A);
    ftk::rand2(b);
    ftk::solve_linear2x2(A, b, x);

    ftk::matrix2x2_vector2_multiplication(A, x, bb);

    REQUIRE(b[0] == Approx(bb[0]).margin(epsilon));
    REQUIRE(b[1] == Approx(bb[1]).margin(epsilon));
  }
}

TEST_CASE("solve_linear3x3") {
  double A[3][3], b[3], x[3], bb[3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand3x3(A);
    ftk::rand3(b);
    ftk::solve_linear3x3(A, b, x);

    ftk::matrix3x3_vector3_multiplication(A, x, bb);

    REQUIRE(b[0] == Approx(bb[0]).margin(epsilon));
    REQUIRE(b[1] == Approx(bb[1]).margin(epsilon));
    REQUIRE(b[2] == Approx(bb[2]).margin(epsilon));
  }
}

TEST_CASE("solve_eigenvalues_symmetric2x2") {
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
      REQUIRE(det == Approx(0.0).margin(epsilon));
    }
  }
}

TEST_CASE("solve_eigenvalues2x2") {
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
      REQUIRE(det == Approx(0.0).margin(epsilon));
    }
  }
}

TEST_CASE("solve_generalized_eigenvalues2x2") {
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
      REQUIRE(det == Approx(0.0).margin(epsilon));
    }
  }
}

TEST_CASE("solve_eigenvectors2x2")
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

      REQUIRE(eig[i] * eigvecs[i][0] == Approx(x[0]).margin(epsilon));
      REQUIRE(eig[i] * eigvecs[i][1] == Approx(x[1]).margin(epsilon));
    }
  }
}

TEST_CASE("solve_eigenvalues_symmetric3x3") {
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
      REQUIRE(det == Approx(0.0).margin(epsilon));
    }
  }
}

TEST_CASE("solve_eigenvectors3x3")
{
  const int nroots = 3; // for now, we use symetric matrix for testing
  double A[3][3], eig[3], eigvecs[3][3], x[3], c[3];
  for (int run = 0; run < nruns; run ++) {
    ftk::rand_symmetric3x3(A);
    ftk::solve_eigenvalues_symmetric3x3(A, eig);
    ftk::solve_eigenvectors3x3(A, nroots, eig, eigvecs);

    for (int i = 0; i < nroots; i ++) {
      ftk::matrix3x3_vector3_multiplication(A, eigvecs[i], x);

      REQUIRE(eig[i] * eigvecs[i][0] == Approx(x[0]).margin(epsilon));
      REQUIRE(eig[i] * eigvecs[i][1] == Approx(x[1]).margin(epsilon));
      REQUIRE(eig[i] * eigvecs[i][2] == Approx(x[2]).margin(epsilon));
    }
  }
}

#if 0 // TODO: check what's wrong here
TEST_CASE(solve_generalized_eigenvalues3x3) {
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
      REQUIRE(0.0, det, epsilon);
    }
  }
}
#endif

int main(int argc, char **argv)
{
  Catch::Session session;
  return session.run(argc, argv);
}
