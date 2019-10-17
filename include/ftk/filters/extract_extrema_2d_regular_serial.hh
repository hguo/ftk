#include "extract_critical_points_2d_regular_serial.hh"
#include <hypermesh/ndarray.hh>
#include <hypermesh/grad.hh>

namespace ftk {

struct extremum_2d_t : public critical_point_2d_t {
  double value;
};

struct extract_extrema_2d_regular_serial : public extract_critical_points_2d_regular_serial {
  extract_extrema_2d_regular_serial() {
    set_symmetric_jacobians(true);
    set_type_filter(ftk::CRITICAL_POINT_2D_MAXIMUM);
  }

  void set_input_scalar_field(const hypermesh::ndarray<double>&);
  void set_input_scalar_field(const double *p, size_t W, size_t H);

  void derive_gradients() {V = hypermesh::gradient2D(scalar);}
  void set_input_gradient_data(const hypermesh::ndarray<double> &grad) {set_input_vector_field(grad);}
  void set_input_gradient_data(const double *p, size_t W, size_t H) {set_input_vector_field(p, W, H);}

  void derive_hessian() {gradV = hypermesh::jacobian2D(V);}
  void set_input_hessian_data(const hypermesh::ndarray<double> &hess) {set_input_jacobian_field(hess);}
  void set_input_hessian_data(const double *p, size_t W, size_t H) {set_input_jacobian_field(p, W, H);}

  void execute();

protected:
  hypermesh::ndarray<double> scalar;
};

/////

void extract_extrema_2d_regular_serial::execute()
{
  if (V.empty()) derive_gradients();
  if (gradV.empty()) derive_hessian();
  if (m.lb() == m.ub()) // unspecified bounds
    m.set_lb_ub({2, 2}, {static_cast<int>(scalar.dim(0)-3), static_cast<int>(scalar.dim(1)-3)});

  fprintf(stderr, "extracting 2D extrema...\n");
  m.element_for(2, [=](hypermesh::regular_simplex_mesh_element e) {
      extremum_2d_t cp;
      if (check_simplex(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        results.push_back(cp);
      }
    }); 
}
  
void extract_extrema_2d_regular_serial::set_input_scalar_field(const hypermesh::ndarray<double>& s)
{
  scalar = s;
}

}
