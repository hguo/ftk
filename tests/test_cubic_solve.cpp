#include <ftk/numeric/cubic_solver.hh>
#include <ftk/numeric/cubic_inequality_solver.hh>
#include <ftk/numeric/cubic_rational_inequality_solver.hh>
#include <iostream>

#if 0
int main(int argc, char **argv)
{
  const double P[4] = {0, 504, -1332, 17388};
  auto roots = ftk::solve_cubic_real_multiplicity(P);
  for (auto kv : roots) 
    fprintf(stderr, "%f, %d\n", kv.first, kv.second);
  return 1;
}
#endif

#if 1
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
  // const float P[4] = {3.f, 1.f, 43.f, 1.f};
  // const float Q[4] = {24.f, -22.f, -4.f, 2.f};

  const double P[4] = {0.000000, 14592.000000, 0.000000, -4435968.000000};
  const double Q[4] = {-32224.000000, -1706080.000000, -7171968.000000, 17743872.000000};
  // float x[3];
  // int n_roots = ftk::solve_cubic_real<float>(P, x);
  // fprintf(stderr, "n_roots=%d, roots={%f, %f, %f}\n", n_roots, x[0], x[1], x[2]);

  // const auto I = ftk::solve_cubic_inequality_real(P);
  // std::cerr << I << std::endl;

  // double x[3], y[3];
  // const int np = ftk::solve_cubic_real(P, x);
  // fprintf(stderr, "np=%d, x={%f, %f, %f}\n", np, x[0], x[1], x[2]);
  // const int nq = ftk::solve_cubic_real(Q, y);
  // fprintf(stderr, "nq=%d, x={%f, %f, %f}\n", np, y[0], y[1], y[2]);

  // const auto I = ftk::solve_cubic_rational_inequality(P, Q);
  const auto I = ftk::solve_cubic_rational_inequality_quantized(P, Q);
  std::cerr << std::get<0>(I) << std::endl;

  return 0;
}

/* expected output (out of date)
x0 = 14.311690 + i*0.000000
x1 = -28.655846 + i*-0.263968
x2 = -28.655846 + i*0.263968
*/

#endif
