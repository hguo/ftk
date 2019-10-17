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

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_results_vtk() const;
#endif

protected:
  hypermesh::ndarray<double> scalar;
  // std::vector<extremum_2d_t> results;
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
        const auto &vertices = e.vertices();
        double values[3];
        for (int i = 0; i < 3; i ++)
          values[i] = scalar(vertices[i][0], vertices[i][1], vertices[i][2]);
        cp.value = lerp_s2(values, cp.mu);
     
        {
          std::lock_guard<std::mutex> guard(mutex);
          results.push_back(cp);
        }
      }
    }); 
}
  
void extract_extrema_2d_regular_serial::set_input_scalar_field(const hypermesh::ndarray<double>& s)
{
  scalar = s;
}

#if FTK_HAVE_VTK
vtkSmartPointer<vtkPolyData> extract_extrema_2d_regular_serial::get_results_vtk() const
{
  fprintf(stderr, "shoot..\n");
  auto polyData = extract_critical_points_2d_regular_serial::get_results_vtk();
#if 0
  vtkSmartPointer<vtkDoubleArray> values = vtkSmartPointer<vtkDoubleArray>::New();
  values->SetNumberOfValues(results.size());
  for (auto i = 0; i < results.size(); i ++) {
    values->SetValue(i, static_cast<double>(results[i].value));
  }
  values->SetName("value");

  polyData->GetPointData()->AddArray(values);
#endif 
  return polyData;
}
#endif

}
