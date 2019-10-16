#include "extract_critical_points_2d_regular_serial.hh"

namespace ftk {

struct extract_extremum_2d_regular_serial : public filter {
  void set_lb_ub(const std::vector<int> &lb, const std::vector<int> &ub);

  void set_input_data(const hypermesh::ndarray<double> &scalar);
  void set_input_data(const double *p, size_t W, size_t H);

  void derive_gradients();
  void set_input_gradient_data(const hypermesh::ndarray<double> &grad);
  void set_input_gradient_data(const double *p, size_t W, size_t H);

  void derive_hessian();
  void set_input_hessian_data(const hypermesh::ndarray<double> &hess);
  void set_input_hessian_data(const double *p, size_t W, size_t H);

protected:
  hypermesh::regular_simplex_mesh m;
  hypermesh::ndarray<double> scalar, grad, hess;

protected:
  bool check_simplex(const hypermesh::regular_simplex_mesh_element& s);
};

/////
void extract_critical_points_2d_regular_serial::set_lb_ub(const std::vector<int> &lb, const std::vector<int> &ub)
{
  m.set_lb_ub(lb, ub);
}

void extract_critical_points_2d_regular_serial::set_input_data(const hypermesh::ndarray<double> &s)
{
  scalar = s;
  set_lb_ub({2, 2}, {s.dim(0)-3, s.dim(1)-3});
}

void extract_critical_points_2d_regular_serial::set_input_data(const double *p, size_t W, size_t H)
{
  scalar = hypermesh::ndarray<double>(p, {W, H});
}

}
