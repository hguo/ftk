#include <iostream>
#include <ftk/numeric/rand.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/linear_solver.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/sign_det.hh>
#include <ftk/ndarray.hh>

int main(int argc, char **argv)
{
#if 0
  long long M[2][2] = {{1, 2}, {1, 2}};

  ftk::print2x2("M", M);
  fprintf(stderr, "det(M)=%lld\n", ftk::det2(M));
  fprintf(stderr, "sign_det(M)=%lld\n", ftk::sign_det2(M));
  fprintf(stderr, "sign_det_delta(M)=%lld\n", ftk::sign_det_delta2(M));
#endif

  // int arr[5] = {1, 5, 4, 3, 2};
  // fprintf(stderr, "nswaps=%d\n", 
  //     ftk::nswaps_bubble_sort<5, int>(arr));

#if 0
  // long long X[3][2] = {{0, 0}, {1, 0}, {0, 1}};
  long long X[3][2] = {{0, 0}, {1, 0}, {-1, 0}};
  int indices[3] = {0, 1, 2};
  fprintf(stderr, "robust_orientation=%lld\n", ftk::robust_orientation2(X, indices));
#endif

#if 0
  // double V[3][2] = {{1.0, 1.0}, {-1.0, 1.0}, {-1.0, -1.0}};
  double V[3][2] = {{1.0, 1.0}, {-1.0, -1.0}, {1.0, -1.0}};
  double mu[3];
  bool b = ftk::inverse_lerp_s2v2(V, mu);
  fprintf(stderr, "b=%d, mu=%f, %f, %f\n", b, mu[0], mu[1], mu[2]);
#endif

  long long V0[3][2] = {{0, 0}, {-1, 1}, {-1, -1}};
  int indices0[3] = {0, 1, 2};
  fprintf(stderr, "b=%d\n", ftk::robust_critical_point_in_simplex2(V0, indices0));
  
  long long V1[3][2] = {{0, 0}, {-1, -1}, {1, -1}};
  int indices1[3] = {0, 2, 3};
  fprintf(stderr, "b=%d\n", ftk::robust_critical_point_in_simplex2(V1, indices1));

  long long V2[3][2] = {{0, 0}, {1, -1}, {1, 0}};
  int indices2[3] = {0, 3, 4};
  fprintf(stderr, "b=%d\n", ftk::robust_critical_point_in_simplex2(V2, indices2));
  
  long long V3[3][2] = {{0, 0}, {1, 0}, {2, -1}};
  int indices3[3] = {0, 4, 5};
  fprintf(stderr, "b=%d\n", ftk::robust_critical_point_in_simplex2(V3, indices3));

  // long long V2[3][2] = {{-1, -1}, {1, -1}, {0, 1}};
  // int indices2[3] = {0, 1, 2};
  // fprintf(stderr, "b=%d\n", ftk::robust_critical_point_in_simplex2(V2, indices2));

  return 0;
}

#if 0
int main(int argc, char **argv)
{
  double X[4][3], V[4][3], J[3][3];
  ftk::rand<double, 4, 3>(X);
  ftk::rand<double, 4, 3>(V);

  ftk::solve_least_square4x3_3(X, V, J);
  ftk::print3x3("J", J);

  return 0;
}
#endif

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
