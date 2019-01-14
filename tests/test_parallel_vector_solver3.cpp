#include <iostream>
#include <vector>
#include <ftk/numeric/parallel_vector_solver.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>

int main(int argc, char **argv)
{
  const double V[3][3] = {
    {0.0000000000, 96.0000000000, 107.0000000000}, 
    {0.0000000000, -147.0000000000, -172.0000000000}, 
    {0.0000000000, -58.0000000000, -147.0000000000}};
  const double W[3][3] = {
    {0.0000000000, -5992.0000000000, 5568.0000000000}, 
    {-3182.0000000000, -9202.0000000000, 27285.5000000000}, 
    {-0.0000000000, -4557.0000000000, 4370.5000000000}};
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
  const auto n = ftk::solve_parallel_vector_barycentric(V, W, lambda, mu);

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
