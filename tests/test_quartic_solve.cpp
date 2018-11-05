#include <ftk/numerics/quartic_solve.hh>
#include <iostream>

int main(int argc, char **argv)
{
  std::complex<float> x[4];
  ftk::quartic_solve<float>(33, 43, 1, 3, x);

  for (int i=0; i<4; i++)
    fprintf(stderr, "x%d = %f + i*%f\n", i, x[i].real(), x[i].imag());

  return 0;
}
