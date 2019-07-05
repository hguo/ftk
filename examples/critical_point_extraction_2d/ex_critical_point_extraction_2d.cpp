#include <mutex>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <hypermesh/ndarray.hh>
#include <hypermesh/regular_simplex_mesh.hh>

hypermesh::ndarray<float> Vf; // vector field
hypermesh::regular_simplex_mesh m(2); // the 2D mesh

std::mutex mutex;

struct critical_point_t {
  float x[2]; // coordinates of the cp
};
 
std::map<hypermesh::regular_simplex_mesh_element, critical_point_t> critical_points;

void check_simplex(const hypermesh::regular_simplex_mesh_element& s)
{
  if (!s.valid()) return; // check if the 2-simplex is valid

  const auto &vertices = s.vertices();
  float X[3][2], v[3][2];

  for (int i = 0; i < 3; i ++) {
    for (int j = 0; j < 2; j ++) {
      v[i][j] = Vf(j, vertices[i][0], vertices[i][1]);
      X[i][j] = vertices[i][j];
    }
  }

  // check intersection
  float mu[3], x[2];
  bool succ = ftk::inverse_lerp_s2v2(v, mu);
  if (!succ) return;

  // check jacobian
  float J[2][2]; // jacobian
  ftk::gradient_2dsimplex2_2(X, v, J);

  std::complex<float> eig[2];
  float delta = ftk::solve_eigenvalues2x2(J, eig);
  
  critical_point_t p;
  ftk::lerp_s2v2(X, mu, x);
  p.x[0] = x[0]; p.x[1] = x[1];
  {
    std::lock_guard<std::mutex> guard(mutex);
    critical_points[s] = p;
  
    std::cerr << s << std::endl;
    fprintf(stderr, "v0={%f, %f}, v1={%f, %f}, v2={%f, %f}, x={%f, %f}\n", 
        v[0][0], v[0][1], v[1][0], v[1][1], v[2][0], v[2][1],
        x[0], x[1]);
    ftk::print2x2("J", J);
    fprintf(stderr, "delta=%f, eig0={%f, %f}, eig1={%f, %f}\n", 
        delta, eig[0].real(), eig[0].imag(), eig[1].real(), eig[1].imag());
  }
}

void extract_critical_points()
{
  fprintf(stderr, "extracting critical points...\n");
  m.set_lb_ub({0, 0}, {static_cast<int>(Vf.dim(1)-1), static_cast<int>(Vf.dim(2)-1)});
  m.element_for(2, check_simplex); // iterate over all 3-simplices
}

void load_vector_field(const std::string &filename)
{
  hypermesh::ndarray<float> U, V;
  U.from_netcdf(filename, "u");
  V.from_netcdf(filename, "v");

  Vf.reshape({2, U.dim(0), U.dim(1)});
  for (int i = 0; i < U.dim(0); i ++) 
    for (int j = 0; j < U.dim(1); j ++) {
      Vf(0, i, j) = U(i, j);
      Vf(1, i, j) = V(i, j);
    }

  fprintf(stderr, "%zu, %zu, %zu\n", Vf.dim(0), Vf.dim(1), Vf.dim(2));
}

int main(int argc, char **argv)
{
  if (argc < 2) return 1;
  load_vector_field(argv[1]);
  extract_critical_points();

  return 0;
}
