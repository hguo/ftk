#include <iostream>
#include <vector>
#include <ftk/numeric/parallel_vector_solver.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>

int main(int argc, char **argv)
{
  const double V[3][3] = {
    {0.0000000000, 0.0000000000, 3.0000000000}, 
    {0.0000000000, 3.0000000000, 0.0000000000}, 
    {0.0000000000, 32.0000000000, 0.0000000000}
  };
  const double W[3][3] = { 
    {0.5000000000, 48.0000000000, 0.0000000000}, 
    {0.5000000000, 48.0000000000, 0.0000000000}, 
    {-608.0000000000, 0.0000000000, 0.0000000000}
  };

#if 0
  const double V[3][3] = {
    {0.0000000000, 3.0000000000, 0.0000000000}, 
    {0.0000000000, 33.0000000000, 0.0000000000}, 
    {0.0000000000, -10.0000000000, 0.0000000000}
  };
  const double W[3][3] = { 
    {-103.5000000000, 49.5000000000, 0.0000000000}, 
    {115.5000000000, 247.5000000000, 0.0000000000}, 
    {400.0000000000, 650.0000000000, 315.0000000000}
  };
#endif

#if 0
  const double V[3][3] = {
    {-65, -54, -59}, 
    {0, -9, 0},
    {0, -124, 0}};
  const double W[3][3] = {
    {6393.5, 1788, -282.5}, 
    {-558, 0, 0}, 
    {-4003, 3348, 4030}};
#endif
#if 0
  const double V[3][3] = {
    {0.0000000000, 0.0000000000, 38.0000000000}, 
    {8.0000000000, 0.0000000000, 0.0000000000}, 
    {3.0000000000, 38.0000000000, 0.0000000000}};
  const double W[3][3] = {
    {608.0000000000, 0.0000000000, 0.0000000000}, 
    {0.0000000000, 0.0000000000, 152.0000000000}, 
    {608.0000000000, 48.0000000000, 152.0000000000}};
#endif

  double lambda[3], mu[3][3];
  const auto n = ftk::solve_parallel_vector_barycentric(V, W, lambda, mu, 1e-10);

  ftk::print3x3("V", V);
  ftk::print3x3("W", W);

  for (int k = 0; k < n; k ++) {
    double v[3], w[3], c[3]; 
    ftk::linear_interpolation3_3(V, mu[k], v);
    ftk::linear_interpolation3_3(W, mu[k], w);
    ftk::cross_product(v, w, c);
    const double norm = ftk::vector_2norm_3(c);
    fprintf(stderr, "lambda=%f, mu={%f, %f, %f}, v={%f, %f, %f}, w={%f, %f, %f}, ||v x w||=%.20f\n", 
        lambda[k], mu[k][0], mu[k][1], mu[k][2], 
        v[0], v[1], v[2], w[0], w[1], w[2], norm); 
  }

  return 0;
}
