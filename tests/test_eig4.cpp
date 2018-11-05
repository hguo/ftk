#include <ftk/numerics/eig.hh>
#include <ftk/numerics/print.hh>
#include <iostream>

int main(int argc, char **argv)
{
  float m[] = {
  0.96170121, 0.25145996, 0.15267866, 0.13361515,
  0.72897646, 0.05956134, 0.59597585, 0.00813389,
  0.82357783, 0.50537611, 0.56218735, 0.77733046,
  0.88185366, 0.51499935, 0.1394362,  0.68746057};
  
  ftk::print4x4("m", m);

  std::complex<float> eig[4], eigvec[4][4];
  ftk::eig4(m, eig, eigvec);

  for (int i=0; i<4; i++) {
    fprintf(stderr, "lambda_%d = %f+%fj\n", i, eig[i].real(), eig[i].imag());
    fprintf(stderr, "vec_%d = (%f+%fj, %f+%fj, %f+%fj, %f+%fj)\n", i,
        eigvec[i][0].real(), eigvec[i][0].imag(),
        eigvec[i][1].real(), eigvec[i][1].imag(),
        eigvec[i][2].real(), eigvec[i][2].imag(),
        eigvec[i][3].real(), eigvec[i][3].imag());
  }

  return 0;
}
