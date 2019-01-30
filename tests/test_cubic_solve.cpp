#include <ftk/numeric/cubic_solver.hh>
#include <ftk/numeric/cubic_inequality_solver.hh>
#include <iostream>

int main(int argc, char **argv)
{
#if 0
  std::complex<float> x[3];
  ftk::solve_cubic<float>(43, 1, 3, x);

  for (int i=0; i<3; i++)
    fprintf(stderr, "x%d = %f + i*%f\n", i, x[i].real(), x[i].imag());

  float P[4] = {3.f, 1.f, 43.f, 1.f};
  ftk::solve_cubic<float>(P, x);
  for (int i=0; i<3; i++)
    fprintf(stderr, "x%d = %f + i*%f\n", i, x[i].real(), x[i].imag());
#endif
  // float P[4] = {24.f, -22.f, -4.f, 2.f};
  float P[4] = {3.f, 1.f, 43.f, 1.f};
  float x[3];
  int n_roots = ftk::solve_cubic_real<float>(P, x);
  fprintf(stderr, "n_roots=%d, roots={%f, %f, %f}\n", n_roots, x[0], x[1], x[2]);

  const auto I = ftk::solve_cubic_inequality_real(P);
  std::cerr << I << std::endl;

  return 0;
}

/* expected output (out of date)
x0 = 14.311690 + i*0.000000
x1 = -28.655846 + i*-0.263968
x2 = -28.655846 + i*0.263968
*/
