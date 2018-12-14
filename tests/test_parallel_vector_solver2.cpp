#include <iostream>
#include <vector>
#include <ftk/numeric/parallel_vector_solver.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>

int main(int argc, char **argv)
{
  for (int i=0; i<10000; i++) {
    double V[4][3], W[4][3];
    double Q[4], P[3][4];

    for (int i = 0; i < 4; i ++) 
      for (int j = 0; j < 3; j ++) {
        V[i][j] = (double)rand() / RAND_MAX - 0.5;
        W[i][j] = (double)rand() / RAND_MAX - 0.5;
      }

    ftk::solve_parallel_vector_tetrahedron(V, W, Q, P);
    fprintf(stderr, "------\n");
    fprintf(stderr, "Q = %f, %f, %f, %f\n", Q[0], Q[1], Q[2], Q[3]);
    for (int i = 0; i < 3; i ++)
      fprintf(stderr, "P[%d] = %f, %f, %f, %f\n", i, P[i][0], P[i][1], P[i][2], P[i][3]);
  }

  return 0;
}
