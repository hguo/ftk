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
#include <ftk/filters/extract_critical_points_2d_regular_serial.hh>
#include <hypermesh/ndarray.hh>
#include <hypermesh/regular_simplex_mesh.hh>

namespace ftk {

struct critical_point_2dt_t {
  double operator[](size_t i) const {if (i > 3) return 0; else return x[i];}
  double x[3];
  double scalar;
  int type = 0;
};

struct track_critical_points_2d_regular_serial : public filter {
  track_critical_points_2d_regular_serial() : m(3) {}
  
  void execute();
  
  void set_input_scalar_field(const double *p, size_t W, size_t H, size_t T);
  void set_input_scalar_field(const hypermesh::ndarray<double>& scalar_) {scalar = scalar_;}

  void set_input_vector_field(const double *p, size_t W, size_t H, size_t T);
  void set_input_vector_field(const hypermesh::ndarray<double>& V_) {V = V_;}

  void set_input_jacobian_field(const double *p, size_t W, size_t H, size_t T); 
  void set_input_jacobian_field(const hypermesh::ndarray<double> &J) {gradV = J;}
  
  void set_lb_ub(const std::vector<int>& lb, const std::vector<int>& ub) {m.set_lb_ub(lb, ub);}

  void get_results();

private:
  hypermesh::ndarray<double> scalar, V, gradV;
  hypermesh::regular_simplex_mesh m;
  
  unsigned int type_filter = 0xffffffff;
  unsigned int jacobian_mode = JACOBIAN_NONE;

  std::map<hypermesh::regular_simplex_mesh_element, critical_point_2dt_t> intersections;
  // std::vector<critical_point_2dt_t> results;

  bool check_simplex(const hypermesh::regular_simplex_mesh_element& s, critical_point_2dt_t& cp);
};


////////////////////
void track_critical_points_2d_regular_serial::execute()
{
  if (!scalar.empty()) {
    if (V.empty()) V = hypermesh::gradient2Dt(scalar);
    if (gradV.empty()) gradV = hypermesh::jacobian2Dt(V);
    jacobian_mode = JACOBIAN_SYMMETRIC;
  }

  if (m.lb() == m.ub()) {
    if (!scalar.empty())
      m.set_lb_ub({2, 2, 0}, {static_cast<int>(V.dim(1)-3), static_cast<int>(V.dim(2)-3), static_cast<int>(V.dim(3)-1)});
    else
      m.set_lb_ub({0, 0, 0}, {static_cast<int>(V.dim(1)-1), static_cast<int>(V.dim(2)-1), static_cast<int>(V.dim(3)-1)});
  }

  fprintf(stderr, "tracking 2D critical points...\n");
  m.element_for(2, [=](hypermesh::regular_simplex_mesh_element e) {
      critical_point_2dt_t cp;
      if (check_simplex(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        intersections[e] = cp;
        // results.push_back(cp);
      }
    }); 
}

bool track_critical_points_2d_regular_serial::check_simplex(
    const hypermesh::regular_simplex_mesh_element& e,
    critical_point_2dt_t& cp)
{
  if (!e.valid()) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(); // obtain the vertices of the simplex
  
  double v[3][2]; // obtain vector values
  for (int i = 0; i < 3; i ++) {
    v[i][0] = V(0, vertices[i][0], vertices[i][1], vertices[i][2]);
    v[i][1] = V(1, vertices[i][0], vertices[i][1], vertices[i][2]);
  }
 
  double mu[3]; // check intersection
  bool succ = inverse_lerp_s2v2(v, mu);
  if (!succ) return false;

  double X[3][3]; // lerp position
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 3; j ++)
      X[i][j] = vertices[i][j];
  lerp_s2v3(X, mu, cp.x);

  return true;

#if 0
  critical_point_2dt_t cp;
  auto eid = f.to_integer();
  ftk::lerp_s2v3(X, mu, I.x);
  I.val = ftk::lerp_s2(value, mu);

  {
    std::lock_guard<std::mutex> guard(mutex);
    intersections[f] = I;
    // fprintf(stderr, "x={%f, %f}, t=%f, val=%f\n", I.x[0], I.x[1], I.x[2], I.val);
  }
#endif
}

}

#endif
