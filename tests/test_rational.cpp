#include <iostream>
#include <ftk/numeric/rational.hh>

int main(int argc, char **argv)
{
  // double Q[4] = {0, 0, -14688, 73008}, 
  //        P[4] = {0, 324, -8586, -4212};

  const double Q[4] = {2000, 12500, 0, 0},
               P[4] = {1000, 0, 0, 0};

  fprintf(stderr, "%f\n", ftk::evaluate_rational<double>(P, Q, 3, -0.16));
  // fprintf(stderr, "%f\n", ftk::evaluate_rational_infinity<double>(P, Q, 3));
  return 0;
}
