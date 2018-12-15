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
    double Q[4], P[4][4];

    for (int i = 0; i < 4; i ++) 
      for (int j = 0; j < 3; j ++) {
        V[i][j] = (double)rand() / RAND_MAX - 0.5;
        W[i][j] = (double)rand() / RAND_MAX - 0.5;
      }

    fprintf(stderr, "====================\n");
    for (int i = 0; i < 4; i ++) {
      double Vp[3][3], Wp[3][3];
      for (int j = 0; j < 3; j ++) {
        int jj = (j >= i) ? (j+1) : j;
        // fprintf(stderr, "face %d, jj=%d\n", i, jj);
        for (int k = 0; k < 3; k ++) {
          Vp[j][k] = V[jj][k];
          Wp[j][k] = W[jj][k];
        }
      }

      double x[3][3];
      const auto n = ftk::solve_parallel_vector_barycentric(Vp, Wp, x);
      if (n > 0) {
        for (int k = 0; k < n; k ++)
          fprintf(stderr, "face %d, x={%f, %f, %f}\n", i, x[k][0], x[k][1], x[k][2]);
      }
    }

#if 1
    ftk::solve_parallel_vector_tetrahedron(V, W, Q, P);
    
    std::complex<double> rootsQ[3], rootsP[4][3];
    ftk::solve_cubic(Q, rootsQ);
    for (int i = 0; i < 4; i ++)
      ftk::solve_cubic(P[i], rootsP[i]);

    fprintf(stderr, "------\n");
    fprintf(stderr, "Q = %f, %f, %f, %f, roots=(%f, %f), (%f, %f), (%f, %f)\n", 
        Q[0], Q[1], Q[2], Q[3], 
        rootsQ[0].real(), rootsQ[0].imag(), 
        rootsQ[1].real(), rootsQ[1].imag(),
        rootsQ[2].real(), rootsQ[2].imag());

    for (int i = 0; i < 4; i ++)
      fprintf(stderr, "P[%d] = %f, %f, %f, %f, roots=(%f, %f), (%f, %f), (%f, %f)\n", i, 
          P[i][0], P[i][1], P[i][2], P[i][3],
          rootsP[i][0].real(), rootsP[i][0].imag(), 
          rootsP[i][1].real(), rootsP[i][1].imag(),
          rootsP[i][2].real(), rootsP[i][2].imag());
    
    // std::vector<std::pair<double, std::array<double, 3> > > pv_points; 
    for (int i = 0; i < 4; i ++) {
      for (int j = 0; j < 3; j ++ ) {
        if (rootsP[i][j].imag() == 0) {
          const double lambda = rootsP[i][j].real();
          const double valQ = ftk::polynomial_evaluate(Q, 3, lambda);
          double x[4]; 
          for (int k = 0; k < 4; k ++) 
            if (k == i) x[k] = 0; // avoid floating error around 0
            else x[k] = ftk::polynomial_evaluate(P[k], 3, lambda) / valQ;
                       // x3 = 1 - x0 - x1 - x2;
          if (x[0] >= 0 && x[0] < 1 && x[1] >= 0 && x[1] < 1 && x[2] >= 0 && x[2] < 1 && x[3] >= 0 && x[3] < 1) {
          // if (1) {
            fprintf(stderr, "lambda=%f, Q(lambda)=%f, x={%f, %f, %f, %f}\n", lambda, valQ, x[0], x[1], x[2], x[3]);
          }
        }
      }
    }
#endif
  }

  return 0;
}
