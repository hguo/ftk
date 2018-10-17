#include <ftk/numerics/cubic_solve.hh>
#include <iostream>

int main(int argc, char **argv)
{
  std::complex<float> x[3];
  ftk::cubic_solve<float>(43, 1, 3, x);

  for (int i=0; i<3; i++)
    fprintf(stderr, "x%d = %f + i*%f\n", i, x[i].real(), x[i].imag());

  return 0;
}

/* expected output:
x0 = 14.311690 + i*0.000000
x1 = -28.655846 + i*-0.263968
x2 = -28.655846 + i*0.263968
*/
