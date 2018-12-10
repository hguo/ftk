#include <ftk/numeric/print.hh>
#include <ftk/numeric/matrix_norm.hh>
#include <ftk/numeric/eigen_solver.hh>
#include <ftk/numeric/cond.hh>
#include <iostream>

int main(int argc, char **argv)
{
  float m[2][2] = {{0.9102779, 0.44108077}, {0.72642273, 0.39278198}};
  ftk::print2x2("m", m);

#if 0
  float m1[2][2] = {{m[0][0], m[1][0]}, {m[0][1], m[1][1]}};
  float m2[2][2];
  ftk::matrix_matrix_multiplication_2x2_2x2(m1, m, m2);
  ftk::print2x2("m^t*m", m2);

  float eig[2];
  ftk::solve_eigenvalues_real_symmetric2(m2, eig);
  fprintf(stderr, "eig=%f, %f\n", eig[0], eig[1]);
#endif

  fprintf(stderr, "1-norm: %f\n", ftk::matrix_1norm_real2x2(m));
  fprintf(stderr, "2-norm: %f\n", ftk::matrix_2norm_real2x2(m));
  fprintf(stderr, "inf-norm: %f\n", ftk::matrix_inf_norm_real2x2(m));
  fprintf(stderr, "fro-norm: %f\n", ftk::matrix_frobenius_norm_real2x2(m));
  fprintf(stderr, "cond: %f\n", ftk::cond_real2x2(m));

  return 0;
}

/* expected output:
m=[[0.9102779031, 0.4410807788], [0.7264227271, 0.3927819729]]
1-norm: 1.636701
2-norm: 1.305495
inf-norm: 1.351359
*/
