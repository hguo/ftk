#include <ftk/numeric/eigen_solver.hh>
#include <ftk/numeric/print.hh>
#include <iostream>

int main(int argc, char **argv)
{
  float A[2][2] = {{0.9102779, 0.44108077}, {0.39278198, 0.95680469}};
  float B[2][2] = {{0.70191435, 0.57763041}, {0.10382177, 0.37101148}};

  ftk::print2x2("A", A);
  ftk::print2x2("B", B);

  std::complex<float> eig[2];
  ftk::solve_generalized_eigenvalues_real2x2(A, B, eig);
  for (int i = 0; i < 2; i ++) {
    fprintf(stderr, "(%f, %f)\n", eig[i].real(), eig[i].imag());
  }
  return 0;
}
