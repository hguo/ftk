#include <iostream>
#include <vector>
#include <ftk/numerics/solve_parallel_vector.hh>
#include <ftk/numerics/cross_product.hh>
#include <ftk/numerics/norm.hh>
#include <ftk/numerics/barycentric_interpolation.hh>

int main(int argc, char **argv)
{
  for (int i=0; i<100; i++) {
    double V[9], W[9], lambda[3] = {0};

    for (int k=0; k<9; k++) {
      V[k] = (double)rand() / RAND_MAX - 0.5;
      W[k] = (double)rand() / RAND_MAX - 0.5;
    }

    auto b = ftk::solve_parallel_vector_barycentric(V, W, lambda);
    if (b) {
      double v[3], w[3], r[3], c[3], cn;
      ftk::barycentric_interpolation3(V, lambda, v);
      ftk::barycentric_interpolation3(W, lambda, w);
      ftk::cross_product(v, w, c);
      cn = ftk::norm2_3(c);

      for (int k=0; k<3; k++) 
        r[k] = v[k] / w[k];

      // if (cn <= 1e-3) { 
      if (1) {
        fprintf(stderr, "[%d] lambda={%f, %f, %f}, v={%f, %f, %f}, w={%f, %f, %f}, ||v x w||=%f\n", i,
            lambda[0], lambda[1], lambda[2],
            v[0], v[1], v[2], w[0], w[1], w[2], cn);
      }
    }
  }

  return 0;
}
