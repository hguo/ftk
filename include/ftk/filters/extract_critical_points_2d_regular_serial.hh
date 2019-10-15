#ifndef _FTK_EXTRACT_CRITICAL_POINT_2D_REGULAR_SERIAL_HH
#define _FTK_EXTRACT_CRITICAL_POINT_2D_REGULAR_SERIAL_HH

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
#include <ftk/filters/filter.hh>

namespace ftk {

struct extract_critical_points_2d_regular_serial : public filter {
  extract_critical_points_2d_regular_serial() : m(2) {}

  void execute();

  void set_input_data(const double *p, size_t W, size_t H) {
    V = hypermesh::ndarray<double>(p, {W, H});
    m.set_lb_ub({0, 0}, {static_cast<int>(W-1), static_cast<int>(H-1)});
  }

  void set_input_data(const hypermesh::ndarray<double> &V_) {
    V = V_;
    m.set_lb_ub({1, 1}, {static_cast<int>(V.dim(1)-2), static_cast<int>(V.dim(2)-2)});
  }
  
  void set_lb_ub(const std::vector<int>& lb, const std::vector<int>& ub) {m.set_lb_ub(lb, ub);}
  
  const std::vector<double>& get_critical_point_coords() const {return coords;}

protected:
  hypermesh::ndarray<double> V;
  hypermesh::regular_simplex_mesh m; // spacetime mesh

  std::vector<double> coords;

  virtual bool check_simplex(const hypermesh::regular_simplex_mesh_element& s);
};

///////

void extract_critical_points_2d_regular_serial::execute()
{
  fprintf(stderr, "extracting 2D critical points...\n");
  m.element_for(2, std::bind(&extract_critical_points_2d_regular_serial::check_simplex, this, std::placeholders::_1)); // iterate over all 3-simplices
}

bool extract_critical_points_2d_regular_serial::check_simplex(const hypermesh::regular_simplex_mesh_element& s)
{
  if (!s.valid()) return false; // check if the 2-simplex is valid

  const auto &vertices = s.vertices();
  double X[3][2], v[3][2];

  for (int i = 0; i < 3; i ++) {
    for (int j = 0; j < 2; j ++) {
      v[i][j] = V(j, vertices[i][0], vertices[i][1]);
      X[i][j] = vertices[i][j];
    }
  }

  // check intersection
  double mu[3], x[2];
  bool succ = ftk::inverse_lerp_s2v2(v, mu);
  if (!succ) return false;

  ftk::lerp_s2v2(X, mu, x);
  {
    std::lock_guard<std::mutex> guard(mutex);
    coords.push_back(x[0]);
    coords.push_back(x[1]);
  }

  return true;
}

}

#endif
