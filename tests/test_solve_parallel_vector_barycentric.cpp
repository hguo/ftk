#include <iostream>
#include <vector>
#include <ftk/numerics/solve_parallel_vector.hh>
#include <ftk/numerics/cross_product.hh>
#include <ftk/numerics/vector_norm.hh>
#include <ftk/numerics/barycentric_interpolation.hh>

int main(int argc, char **argv)
{
  for (int i=0; i<10000; i++) {
    double V[9], W[9], lambda[9] = {0};

    for (int k=0; k<9; k++) {
      V[k] = (double)rand() / RAND_MAX - 0.5;
      W[k] = (double)rand() / RAND_MAX - 0.5;
    }

    const auto n = ftk::solve_parallel_vector_barycentric(V, W, lambda);
    for (int j = 0; j < n; j ++) {
      double v[3], w[3], r[3], c[3], cn;
      ftk::barycentric_interpolation3(V, lambda+j*3, v);
      ftk::barycentric_interpolation3(W, lambda+j*3, w);
      ftk::cross_product(v, w, c);
      cn = ftk::vector_norm2_3(c);

      for (int k=0; k<3; k++) 
        r[k] = v[k] / w[k];

      // if (cn <= 1e-3) { 
      if (1) {
        fprintf(stderr, "[%d] lambda={%f, %f, %f}, v={%f, %f, %f}, w={%f, %f, %f}, ||v x w||=%.20f\n", i,
            lambda[j*3+0], lambda[j*3+1], lambda[j*3+2],
            v[0], v[1], v[2], w[0], w[1], w[2], cn);
      }
    }
  }

  return 0;
}
