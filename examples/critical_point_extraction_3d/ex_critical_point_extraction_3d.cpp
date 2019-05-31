#include <mutex>

#if HAVE_NETCDF
#include <netcdf.h>
#endif

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

const int DW = 128, DH = 128, DD = 128;// the dimensionality of the data is DW*DH

hypermesh::ndarray<float> scalar, grad, hess;
hypermesh::regular_simplex_mesh m(3); // the 3D spatial mesh

std::mutex mutex;

struct critical_point_t {
  float x[4]; // the spacetime coordinates of the trajectory
};
 
std::map<hypermesh::regular_simplex_mesh_element, critical_point_t> critical_points;

// the output trajectories
std::vector<std::vector<std::vector<float>>> trajectories;

void derive_gradients()
{
  fprintf(stderr, "deriving gradients...\n");
  grad.reshape({3, DW, DH, DD});
  for (int k = 1; k < DD-1; k ++) {
    for (int j = 1; j < DH-1; j ++) {
      for (int i = 1; i < DW-1; i ++) {
        grad(0, i, j, k) = 0.5 * (scalar(i+1, j, k) - scalar(i-1, j, k));
        grad(1, i, j, k) = 0.5 * (scalar(i, j+1, k) - scalar(i, j-1, k));
        grad(2, i, j, k) = 0.5 * (scalar(i, j, k+1) - scalar(i, j, k-1));
      }
    }
  }
}

void derive_hessians()
{
  fprintf(stderr, "deriving hessians...\n");
  hess.reshape({3, 3, DW, DH, DD});

  for (int k = 0; k < DD-2; k ++) {
    for (int j = 2; j < DH-2; j ++) {
      for (int i = 2; i < DW-2; i ++) {
        const float H00 = hess(0, 0, i, j, k) = // ddf/dx2
          0.5 * (grad(0, i+1, j, k) - grad(0, i-1, j, k));
        const float H01 = hess(0, 1, i, j, k) = // ddf/dxdy
          0.5 * (grad(0, i, j+1, k) - grad(0, i, j-1, k));
        const float H02 = hess(0, 2, i, j, k) = // ddf/dxdz
          0.5 * (grad(0, i, j, k+1) - grad(0, i, j, k-1));

        const float H10 = hess(1, 0, i, j, k) = // ddf/dydx
          0.5 * (grad(1, i+1, j, k) - grad(1, i-1, j, k));
        const float H11 = hess(1, 1, i, j, k) = // ddf/dy2
          0.5 * (grad(1, i, j+1, k) - grad(1, i, j-1, k));
        const float H12 = hess(1, 2, i, j, k) = // ddf/dydz
          0.5 * (grad(1, i, j, k+1) - grad(1, i, j, k-1));

        const float H20 = hess(2, 0, i, j, k) = // ddf/dydx
          0.5 * (grad(2, i+1, j, k) - grad(2, i-1, j, k));
        const float H21 = hess(2, 1, i, j, k) = // ddf/dy2
          0.5 * (grad(2, i, j+1, k) - grad(2, i, j-1, k));
        const float H22 = hess(2, 2, i, j, k) = // ddf/dydz
          0.5 * (grad(2, i, j, k+1) - grad(2, i, j, k-1));
      }
    }
  }
}

void check_simplex(const hypermesh::regular_simplex_mesh_element& s)
{
  if (!s.valid()) return; // check if the 3-simplex is valid
  // fprintf(stderr, "%zu\n", s.to_integer());

  const auto &vertices = s.vertices();
  float X[4][4], g[4][3], value[4];

  for (int i = 0; i < 4; i ++) {
    for (int j = 0; j < 3; j ++)
      g[i][j] = grad(j, vertices[i][0], vertices[i][1], vertices[i][2]); // , vertices[i][3]);
    for (int j = 0; j < 4; j ++)
      X[i][j] = vertices[i][j];
    value[i] = scalar(vertices[i][0], vertices[i][1], vertices[i][2]); // , vertices[i][3]);
  }

  // check intersection
  float mu[4], x[3];
  bool succ = ftk::inverse_lerp_s3v3(g, mu);
  if (!succ) return;

  // check hessian
  float H[4][3][3], h[3][3];
  for (int i = 0; i < 4; i ++)
    for (int j = 0; j < 3; j ++)
      for (int k = 0; k < 3; k ++)
        H[i][j][k] = hess(j, k, vertices[i][0], vertices[i][1], vertices[i][2]); // , vertices[i][3]);
  ftk::lerp_s3m3x3(H, mu, h);

  float eig[3];
  ftk::solve_eigenvalues_symmetric3x3(h, eig);
  // fprintf(stderr, "eig=%f, %f, %f\n", eig[0], eig[1], eig[2]);

  if (eig[0] < 0 && eig[1] < 0 && eig[2] < 0) { // local maxima
    // dump results
    float val = ftk::lerp_s3(value, mu);
    ftk::lerp_s3v4(X, mu, x);
   
    critical_point_t p;
    p.x[0] = x[0]; p.x[1] = x[1]; p.x[2] = x[2]; // p.x[3] = x[3];
    {
      std::lock_guard<std::mutex> guard(mutex);
      critical_points[s] = p;
    
      std::cerr << s << std::endl;
      fprintf(stderr, "x={%f, %f, %f}\n", x[0], x[1], x[2]); // , x[3]);
    }
  }
}

void extract_critical_points()
{
  fprintf(stderr, "extracting critical points...\n");
  m.set_lb_ub({2, 2, 2}, {DW-3, DH-3, DD-3}); // set the lower and upper bounds of the mesh
  m.element_for(3, check_simplex); // iterate over all 3-simplices
}

int main(int argc, char **argv)
{
  size_t starts[4] = {0, 0, 0, 0}, 
         sizes[4]  = {1, size_t(DD), size_t(DH), size_t(DW)};

  scalar.from_netcdf(argv[1], "vort", starts, sizes);
  scalar.reshape({DW, DH, DD});

  derive_gradients();
  derive_hessians();
  
  extract_critical_points();

  return 0;
}
