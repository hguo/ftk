#include <iostream>
#include <ftk/numeric/rand.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/linear_solver.hh>
#include <ftk/ndarray.hh>

int main(int argc, char **argv)
{
  double X[4][3], V[4][3], J[3][3];
  ftk::rand<double, 4, 3>(X);
  ftk::rand<double, 4, 3>(V);

  ftk::solve_least_square4x3_3(X, V, J);
  ftk::print3x3("J", J);

  return 0;
}

#if 0
int main(int argc, char **argv)
{
  const int m = 7, n = 8; 
  ftk::ndarray<int> array({m, n}); // n rows and m columns
  for (int i = 0; i < array.nelem(); i ++)
    array[i] = i;
  
  for (int j = 0; j < n; j ++) {
    for (int i = 0; i < m; i ++) 
      fprintf(stderr, "%d, ", array(i, j));
    fprintf(stderr, "\n");
  }

  const int p = 3, q = 2;
  auto subarray = array.slice({3, 2}, {p, q});
  
  for (int j = 0; j < q; j ++) {
    for (int i = 0; i < p; i ++) 
      fprintf(stderr, "%d, ", subarray(i, j));
    fprintf(stderr, "\n");
  }

  return 0;
}
#endif
