#include <ftk/numeric/eigen_solver.hh>
#include <ftk/numeric/print.hh>
#include <iostream>

int main(int argc, char **argv)
{
  float A[3][3] = {{0.9102779, 0.44108077, 0.72642273}, {0.39278198, 0.95680469, 0.02683596}, {0.05335823, 0.86960914, 0.43971526}};
  float B[3][3] = {{0.70191435, 0.57763041, 0.59147791}, {0.10382177, 0.37101148, 0.30612991}, {0.86739255, 0.41015195, 0.45043249}};

  ftk::print3x3("A", A);
  ftk::print3x3("B", B);

  std::complex<float> eig[3];
  ftk::solve_generalized_eigenvalues_real3x3(A, B, eig);
  for (int i = 0; i < 3; i ++) {
    fprintf(stderr, "(%f, %f)\n", eig[i].real(), eig[i].imag());
  }

#if 0
  float m[] = {0.9102779, 0.44108077, 0.72642273, 0.39278198, 0.95680469, 0.02683596, 0.05335823, 0.86960914, 0.43971526};
  std::complex<float> eig[3], eigvec[3][3];

  ftk::eig3(m, eig, eigvec);

  for (int i=0; i<3; i++) {
    fprintf(stderr, "lambda_%d = %f+%fj\n", i, eig[i].real(), eig[i].imag());
    fprintf(stderr, "vec_%d = (%f+%fj, %f+%fj, %f+%fj)\n", i,
        eigvec[i][0].real(), eigvec[i][0].imag(),
        eigvec[i][1].real(), eigvec[i][1].imag(),
        eigvec[i][2].real(), eigvec[i][2].imag());
  }
 
  fprintf(stderr, "=====\n");

  float m1[3][3] = {{0.9102779, 0.44108077, 0.72642273}, {0.39278198, 0.95680469, 0.02683596}, {0.05335823, 0.86960914, 0.43971526}};
  ftk::print3x3("m1", m1);
  ftk::eig3(m1, eig, eigvec);

  for (int i=0; i<3; i++) {
    fprintf(stderr, "lambda_%d = %f+%fj\n", i, eig[i].real(), eig[i].imag());
    fprintf(stderr, "vec_%d = (%f+%fj, %f+%fj, %f+%fj)\n", i,
        eigvec[i][0].real(), eigvec[i][0].imag(),
        eigvec[i][1].real(), eigvec[i][1].imag(),
        eigvec[i][2].real(), eigvec[i][2].imag());
    // fprintf(stderr, "nrom=%f\n", ftk::vector_2norm_3(eigvec[2]).real());
  }

#endif
  return 0;
}

/* expected output: 
lambda_0 = 1.586087+0.000000j
vec_0 = (1.860765+0.000000j, 1.204087+0.000000j, 1.000000+0.000000j)
lambda_1 = 0.360355+0.428522j
vec_1 = (-0.581230+-0.892065j, -0.055596+0.547512j, 1.000000+0.000000j)
lambda_2 = 0.360355+-0.428522j
vec_2 = (-0.581230+0.892065j, -0.055596+-0.547512j, 1.000000+0.000000j)
*/
