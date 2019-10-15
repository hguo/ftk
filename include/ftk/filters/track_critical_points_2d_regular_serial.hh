#ifndef _FTK_TRACK_CRITICAL_POINTS_2D_REGULAR_SERIAL_HH
#define _FTK_TRACK_CRITICAL_POINTS_2D_REGULAR_SERIAL_HH

#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
// #include <ftk/algorithms/cca.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <hypermesh/ndarray.hh>
#include <hypermesh/regular_simplex_mesh.hh>

namespace ftk {

struct track_critical_points_2d_regular_serial {
  track_critical_points_2d_regular_serial() : m(3) {}
  
  void execute();

  void set_input_data(const double *p, size_t W, size_t H, size_t T) {
    V = hypermesh::ndarray<double>(p, {W, H, T});
    m.set_lb_ub({0, 0, 0}, {static_cast<int>(W), static_cast<int>(H), static_cast<int>(T)});
  }

  void set_input_data(const hypermesh::ndarray<double> &V_) {
    V = V_;
    m.set_lb_ub({0, 0, 0}, {static_cast<int>(V.dim(1)-2), static_cast<int>(Vf.dim(2)-2), static_cast<int>(Vf.dim(3)-2)});
  }
  
  void set_lb_ub(const std::vector<int>& lb, const std::vector<int>& ub) {m.set_lb_ub(lb, ub);}

  void get_results();

private:
  hypermesh::ndarray<double> V;
  hypermesh::regular_simplex_mesh m;
};


////////////////////
void check_simplex(const hypermesh::regular_simplex_mesh_element& f)
{
  if (!f.valid()) return; // check if the 2-simplex is valid
  const auto &vertices = f.vertices(); // obtain the vertices of the simplex
  double X[3][2], v[3][2];

  for (int i = 0; i < 3; i ++) {
    v[i][0] = V(0, vertices[i][0], vertices[i][1], vertices[i][2]);
    v[i][1] = V(1, vertices[i][0], vertices[i][1], vertices[i][2]);
  }
 
  float mu[3];
  bool succ = ftk::inverse_lerp_s2v2(v, mu);
  if (!succ) return;

  for (int i = 0; i < vertices.size(); i ++)
    for (int j = 0; j < 3; j ++)
      X[i][j] = vertices[i][j];

  intersection_t I;
  I.eid = f.to_integer();
  ftk::lerp_s2v3(X, mu, I.x);
  I.val = ftk::lerp_s2(value, mu);

  {
    std::lock_guard<std::mutex> guard(mutex);
    intersections[f] = I;
    // fprintf(stderr, "x={%f, %f}, t=%f, val=%f\n", I.x[0], I.x[1], I.x[2], I.val);
  }
}

}

#endif
