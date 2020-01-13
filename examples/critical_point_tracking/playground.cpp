#include <iostream>
#include <ftk/numeric/rand.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/linear_solver.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/sign_det.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/adjugate.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/rational_number.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>

#if 0
int main(int argc, char **argv)
{
  // fprintf(stderr, "%d\n", ftk::gcd(100, -10));
#if 1 // test adjugate 3x3
  long long A[3][3] = {
    {-3, 2, -5}, 
    {-1, 0, -2}, 
    {3, -4, 1}
  };
  long long B[3][3];
  ftk::adjugate3(A, B);
  ftk::print3x3("B", B);
#endif

#if 0
  ftk::rational<long long> a(3, 4), b(5, 9);
  std::cerr << a - b << std::endl;
#endif

  return 0;
}
#endif

#if 1 // cp2d extraction
const int DW = 63, DH = 63;
const long long factor = 1000;

ftk::ndarray<double> Vf; // vector field
ftk::regular_simplex_mesh m(2); // the 2D mesh

std::mutex mutex;

struct critical_point_t {
  double x[2]; // coordinates of the cp
};
 
std::map<ftk::regular_simplex_mesh_element, critical_point_t> critical_points;

void check_simplex0(const ftk::regular_simplex_mesh_element& s)
{
  if (!s.valid(m)) return; 
  const auto &vertices = s.vertices(m);

  double v[2];
  long long iv[2];
  for (int j = 0; j < 2; j ++) {
    v[j] = Vf(j, vertices[0][1], vertices[0][1]);
    iv[j] = v[j] * factor;
  }

  if (iv[0] == 0 && iv[1] == 0) {
    fprintf(stderr, "found cp on 0-element\n");
  }
}

void check_simplex1(const ftk::regular_simplex_mesh_element& s)
{
  if (!s.valid(m)) return; 
  const auto &vertices = s.vertices(m);

  double V[2][2], X[2][2];
  long long iV[2][2];
  for (int i = 0; i < 2; i ++) {
    for (int j = 0; j < 2; j ++) {
      V[i][j] = Vf(j, vertices[i][0], vertices[i][1]);
      iV[i][j] = V[i][j] * factor;
      X[i][j] = vertices[i][j];
    }
  }

  bool succ = ftk::integer_inverse_lerp_s1v2(iV);
  if (succ) {
#if 0
    fprintf(stderr, "vertices=%d, %d, %d, %d\n",
        vertices[0][0], vertices[0][1],
        vertices[1][0], vertices[1][1]);
    ftk::print2x2("V", V);
    ftk::print2x2("iV", iV);
    fprintf(stderr, "det=%lld\n", ftk::det2(iV));
#endif

    double mu;
    ftk::inverse_lerp_s1v2(V, mu);
    double x[2];
    ftk::lerp_s1v2(X, mu, x);
    fprintf(stderr, "found cp on 1-element, mu=%f, x={%f, %f}\n", 
        mu, x[0], x[1]);
  }
}

void check_simplex(const ftk::regular_simplex_mesh_element& s)
{
  if (!s.valid(m)) return; // check if the 2-simplex is valid

  const auto &vertices = s.vertices(m);
  double X[3][2], v[3][2];
  long long intv[3][2];

  for (int i = 0; i < 3; i ++) {
    for (int j = 0; j < 2; j ++) {
      v[i][j] = Vf(j, vertices[i][0], vertices[i][1]);
      intv[i][j] = factor * v[i][j];
      X[i][j] = vertices[i][j];
    }
  }
#if 0 // quantization based detection
  const long long M[3][3] = {
    {intv[0][0], intv[1][0], intv[2][0]},
    {intv[1][1], intv[1][1], intv[2][1]},
    {factor, factor, factor}
  };
  long long adjM[3][3];
  long long detM = ftk::det3(M);
  ftk::adjugate3(M, adjM);
  long long b[3] = {adjM[0][2]*factor, adjM[1][2]*factor, adjM[2][2]*factor};
  // long long b[3] = {adjM[0][2], adjM[1][2], adjM[2][2]};

  int sign_detM = ftk::sign(detM);
  if (sign_detM < 0) {
    detM *= sign_detM;
    for (int k = 0; k < 3; k ++)
      b[k] *= sign_detM;
  }
  bool succ = (b[0] > 0 && b[0] < detM && b[1] > 0 && b[1] < detM && b[2] > 0 && b[2] < detM);
  // bool succ = (b[0] >= 0 && b[0] <= detM && b[1] >= 0 && b[1] <= detM && b[2] >= 0 && b[2] <= detM);

  if (!succ) return;
#if 0
  fprintf(stderr, "detM=%lld, b=%lld, %lld, %lld, m=%f, %f, %f\n", detM, b[0], b[1], b[2], 
      (double)b[0]/detM, (double)b[1]/detM, (double)b[2]/detM);
#endif
  // fprintf(stderr, "det=%lld, b=%lld, %lld, %lld\n", 
  //     detM, b[0], b[1], b[2]);

  // pre-check with robust cp detection
  int indices[3];
  for (int i = 0; i < 3; i ++) {
    std::vector<int> my_vertices = {vertices[i][0], vertices[i][1], vertices[i][2]};
    indices[i] = m.get_lattice().to_integer(vertices[i]);
  }

  // const int sign = ftk::positive2(X, indices);
  // bool succ0 = ftk::robust_critical_point_in_simplex2(intv, indices, sign);
  // if (!succ0) return;
#endif
  // check intersection
  double mu[3], x[2];
  bool succ1 = ftk::inverse_lerp_s2v2(v, mu);
  // fprintf(stderr, "indices=%d, %d, %d, mu=%f, %f, %f, ", 
  //     indices[0], indices[1], indices[2], mu[0], mu[1], mu[2]);
  // fprintf(stderr, "sign=%d, indices=%d, %d, %d, ", 
  //     sign, indices[0], indices[1], indices[2]);
  if (!succ1) return;

  // check jacobian
  double J[2][2]; // jacobian
  // ftk::gradient_2dsimplex2_2(X, v, J);
  ftk::jacobian_2dsimplex2(X, v, J);

  std::complex<double> eig[2];
  double delta = ftk::solve_eigenvalues2x2(J, eig);
  
  critical_point_t p;
  ftk::lerp_s2v2(X, mu, x);
  p.x[0] = x[0]; p.x[1] = x[1];
  {
    std::lock_guard<std::mutex> guard(mutex);
    critical_points[s] = p;
  
    // std::cerr << s << std::endl;
    fprintf(stderr, "v0={%f, %f}, v1={%f, %f}, v2={%f, %f}, mu={%f, %f, %f}, x={%f, %f}\n", 
        v[0][0], v[0][1], v[1][0], v[1][1], v[2][0], v[2][1],
        mu[0], mu[1], mu[2], x[0], x[1]);
    // fprintf(stderr, "v0={%lld, %lld}, v1={%lld, %lld}, v2={%lld, %lld}, x={%f, %f}\n", 
    //     intv[0][0], intv[0][1], intv[1][0], intv[1][1], intv[2][0], intv[2][1],
    //     x[0], x[1]);
#if 0
    ftk::print2x2("J", J);
    fprintf(stderr, "delta=%f, eig0={%f, %f}, eig1={%f, %f}\n", 
        delta, eig[0].real(), eig[0].imag(), eig[1].real(), eig[1].imag());

    if (delta >= 0) { // two real roots
      if (eig[0].real() * eig[1].real() < 0) {
        fprintf(stderr, "SADDLE\n");
      } else if (eig[0].real() * eig[1].real() > 0) {
        fprintf(stderr, "ATTRACTING/REPELLING\n");
      }
    } else { // two conjugate roots
      if (eig[0].real() < 0) {
        fprintf(stderr, "ATTRACTING FOCUS\n");
      } else if (eig[0].real() > 0) {
        fprintf(stderr, "REPELLING FOCUS\n");
      } else 
        fprintf(stderr, "CENTER\n");
    }
#endif
  }
}

void extract_critical_points()
{
  fprintf(stderr, "extracting critical points...\n");
  m.set_lb_ub({2, 2}, {DW-3, DH-3}); // {static_cast<int>(Vf.dim(1)-1), static_cast<int>(Vf.dim(2)-1)});
  // m.element_for(2, check_simplex, 1); // iterate over all 3-simplices
  m.element_for(0, check_simplex0, 1); // iterate over all 3-simplices
  m.element_for(1, check_simplex1, 1); // iterate over all 3-simplices
  m.element_for(2, check_simplex, 1); // iterate over all 3-simplices
}

int main(int argc, char **argv)
{
  auto scalar = ftk::synthetic_woven_2D<double>(DW, DH, 1e-3); // , 0.0);
  Vf = gradient2D(scalar);
  extract_critical_points();

  return 0;
}
#endif // cp2d extraction


#if 0 // misc test
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

  // two test case from woven
#if 0
  // long long V0[3][2] = {{17751, 7354}, {0, 7406}, {0, -10396}};
  long long V0[3][2] = {{1, 1}, {0, 1}, {0, -1}};
  // int indices0[3] = {7341, 7463, 7464}; 
  int indices0[3] = {0, 1, 2}; 
  fprintf(stderr, "----b=%d\n", ftk::robust_critical_point_in_simplex2(V0, indices0));
#endif
 
#if 0
  // long long V1[3][2] = {{0, 7406}, {-17729, -10323}, {0, -10396}}; 
  // int indices1[3] = {7463, 7586, 7464};
  long long V1[3][2] = {{0, -1}, {0, 1}, {-1, -1}}; 
  // long long V1[3][2] = {{0, 7406}, {0, -10396}, {-17729, -10323}}; 
  // int indices1[3] = {7463, 7464, 7586};
  int indices1[3] = {2, 1, 3}; 
  fprintf(stderr, "----b=%d\n", ftk::robust_critical_point_in_simplex2(V1, indices1));
#endif

#if 1
  // point in simplex test
  long long X0[3][2] = {{1, 1}, {0, 1}, {0, -1}};
  int indices0[3] = {0, 1, 2};
  long long X1[3][2] = {{0, -1}, {0, 1}, {-1, -1}}; 
  int indices1[3] = {2, 1, 3}; 
  long long zero[2] = {0, 0};

  fprintf(stderr, "--b=%d\n", ftk::robust_point_in_simplex2(X0, indices0, zero));
  fprintf(stderr, "--b=%d\n", ftk::robust_point_in_simplex2(X1, indices1, zero));
  // fprintf(stderr, "--b=%d\n", ftk::robust_point_in_polygon2(zero, 3, indices0, X0));
  // fprintf(stderr, "--b=%d\n", ftk::robust_point_in_polygon2(zero, 3, indices1, X1));
#endif

#if 0
  long long V0[3][2] = {{0, 0}, {-1, 1}, {-1, -1}};
  int indices0[3] = {0, 1, 2};
  fprintf(stderr, "---b=%d\n", ftk::robust_critical_point_in_simplex2(V0, indices0));
  
  long long V1[3][2] = {{0, 0}, {-1, -1}, {1, -1}};
  int indices1[3] = {0, 2, 3};
  fprintf(stderr, "---b=%d\n", ftk::robust_critical_point_in_simplex2(V1, indices1));

  long long V2[3][2] = {{0, 0}, {1, -1}, {1, 0}};
  int indices2[3] = {0, 3, 4};
  fprintf(stderr, "---b=%d\n", ftk::robust_critical_point_in_simplex2(V2, indices2));
  
  long long V3[3][2] = {{0, 0}, {1, 0}, {2, -1}};
  int indices3[3] = {0, 4, 5};
  fprintf(stderr, "---b=%d\n", ftk::robust_critical_point_in_simplex2(V3, indices3));
#endif

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
#endif
