#include <ftk/numerics/eig.hh>
#include <iostream>

int main(int argc, char **argv)
{
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
