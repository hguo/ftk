#include <ftk/numerics/eig.hh>
#include <iostream>

int main(int argc, char **argv)
{
  float m[] = {0.9102779, 0.44108077, 0.72642273, 0.39278198, 0.95680469, 0.02683596, 0.05335823, 0.86960914, 0.43971526};
  std::complex<float> eig[3], eigvec[3][3];

  ftk::eig3(m, eig, eigvec);

  for (int i=0; i<3; i++)
    fprintf(stderr, "lambda_%d = %f + i*%f\n", i, eig[i].real(), eig[i].imag());

  return 0;
}
