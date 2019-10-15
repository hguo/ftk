#include "extract_critical_points_2d_regular_serial.hh"

namespace ftk {

struct extract_extremum_2d_regular_serial : public filter {
  void set_input_scalar_data(const hypermesh::ndarray<double> &scalar);
  void set_input_scalar_data(const double *p, size_t W, size_t H);

  void derive_gradients();
  void set_input_gradient_data(const hypermesh::ndarray<double> &grad);
  void set_input_gradient_data(const double *p, size_t W, size_t H);

  void derive_hessian();
  void set_input_hessian_data(const hypermesh::ndarray<double> &hess);
  void set_input_hessian_data(const double *p, size_t W, size_t H);
};

}
