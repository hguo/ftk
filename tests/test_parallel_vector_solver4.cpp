#include <iostream>
#include <vector>
#include <ftk/numeric/parallel_vector_solver.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>

int main(int argc, char **argv)
{
  double V[4][3], W[4][3];
#if 0
  V[0][0] = 0;  V[0][1] = 0;  V[0][2] = 38;
  V[1][0] = 8;  V[1][1] = 0;  V[1][2] = 0;
  V[2][0] = 3;  V[2][1] = 38;  V[2][2] = 0;
  V[3][0] = 32;  V[3][1] = 0;  V[3][2] = -8;

  W[0][0] = 608;  W[0][1] = 0;  W[0][2] = 0;
  W[1][0] = 0;  W[1][1] = 0;  W[1][2] = 152;
  W[2][0] = 608;  W[2][1] = 48;  W[2][2] = 152; 
  W[3][0] = -1216;  W[3][1] = 48;  W[3][2] = 0;
#endif
  
  V[0][0] = 0;  V[0][1] = 8;  V[0][2] = 0;
  V[1][0] = 0;  V[1][1] = 0;  V[1][2] = 0;
  V[2][0] = 0;  V[2][1] = 32;  V[2][2] = 0;
  V[3][0] = 3;  V[3][1] = 38;  V[3][2] = 0;

  W[0][0] = 0;  W[0][1] = 0;  W[0][2] = 152;
  W[1][0] = 0;  W[1][1] = 0;  W[1][2] = 0;
  W[2][0] = -608;  W[2][1] = 0;  W[2][2] = 0; 
  W[3][0] = 608;  W[3][1] = 48;  W[3][2] = 152;

  // auto I = ftk::solve_parallel_vector_tetrahedron_inequalities(V, W);
  auto I = ftk::solve_parallel_vector_tetrahedron_inequalities_quantized(V, W);
  std::cerr << I << std::endl;
  
  return 0;
}
