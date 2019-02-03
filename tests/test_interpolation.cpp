#include <ftk/numeric/linear_interpolation.hh>
#include <iostream>

int main(int argc, char **argv)
{
  const double X[4][3] = {
    {95, 115, 9},
    {95, 114, 10}, 
    {96, 115, 10}, 
    {95, 115, 10}
  };
  const double mu[4] = {0, 0, 0, 1};
  double x[3]; 
  
  ftk::linear_interpolation4_3(X, mu, x);
  fprintf(stderr, "x={%f, %f, %f}\n", x[0], x[1], x[2]);

  return 0;
}
