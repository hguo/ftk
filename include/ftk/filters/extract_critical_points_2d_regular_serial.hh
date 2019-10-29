#ifndef _FTK_EXTRACT_CRITICAL_POINT_2D_REGULAR_SERIAL_HH
#define _FTK_EXTRACT_CRITICAL_POINT_2D_REGULAR_SERIAL_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/critical_point.hh>
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
#include <ftk/ndarray.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/hypermesh/regular_simplex_mesh.hh>
#include <ftk/filters/filter.hh>

#if FTK_HAVE_VTK
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>
#endif

namespace ftk {

struct critical_point_2d_t {
  double operator[](size_t i) const {if (i > 2) return 0; else return x[i];}

  double x[2];
  double scalar; // used only if the scalar field is available
  int type = 0;
};

struct extract_critical_points_2d_regular_serial : public filter {
  extract_critical_points_2d_regular_serial() : m(2) {}

  void execute();

  void set_input_scalar_field(const double *p, size_t W, size_t H);
  void set_input_scalar_field(const ndarray<double>&);

  void set_input_vector_field(const double *p, size_t W, size_t H);
  void set_input_vector_field(const ndarray<double>&);

  void set_input_jacobian_field(const double *p, size_t W, size_t H); // must be the same dimension as the input data
  void set_input_jacobian_field(const ndarray<double> &J) {gradV = J;}

  void set_lb_ub(const std::vector<int>& lb, const std::vector<int>& ub) {m.set_lb_ub(lb, ub);}
  void set_type_filter(unsigned int mask = 0xffffffff) {type_filter = mask;}

  virtual const std::vector<critical_point_2d_t>& get_results() const {return results;}

#if FTK_HAVE_VTK
  void set_input_scalar_field(vtkSmartPointer<vtkImageData>);
  vtkSmartPointer<vtkPolyData> get_results_vtk() const;
#endif

protected:
  ndarray<double> scalar, V, gradV;
  regular_simplex_mesh m; // spacetime mesh
  
  unsigned int type_filter = 0xffffffff;
  bool has_jacobian = false;
  bool symmetric_jacobian = false;

  std::vector<critical_point_2d_t> results;

  bool check_simplex(const regular_simplex_mesh_element& s, critical_point_2d_t& cp);
};

///////

void extract_critical_points_2d_regular_serial::set_input_scalar_field(const ndarray<double>& s)
{
  scalar = s;
}
  
void extract_critical_points_2d_regular_serial::set_input_vector_field(const double *p, size_t W, size_t H)
{
  V = ndarray<double>(p, {2, W, H});
}

void extract_critical_points_2d_regular_serial::set_input_vector_field(const ndarray<double> &V_) 
{
  V = V_;
}

void extract_critical_points_2d_regular_serial::set_input_jacobian_field(const double *p, size_t W, size_t H)
{
  gradV = ndarray<double>(p, {2, 2, W, H});
}

void extract_critical_points_2d_regular_serial::execute()
{
  if (!scalar.empty()) {
    if (V.empty()) V = gradient2D(scalar);
    if (gradV.empty()) gradV = jacobian2D(V);
    has_jacobian =  true;
    symmetric_jacobian = true;
  }

  if (m.lb() == m.ub()) { // unspecified bounds
    if (!scalar.empty()) 
      m.set_lb_ub({2, 2}, {static_cast<int>(V.dim(1)-3), static_cast<int>(V.dim(2)-3)});
    else
      m.set_lb_ub({0, 0}, {static_cast<int>(V.dim(1)-1), static_cast<int>(V.dim(2)-1)});
  }

  fprintf(stderr, "extracting 2D critical points...\n");
  m.element_for(2, [=](regular_simplex_mesh_element e) {
      critical_point_2d_t cp;
      if (check_simplex(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        results.push_back(cp);
      }
    }); 
}

bool extract_critical_points_2d_regular_serial::check_simplex(
    const regular_simplex_mesh_element& s, 
    critical_point_2d_t &cp)
{
  if (!s.valid()) return false; // check if the 2-simplex is valid
  const auto &vertices = s.vertices();

  double v[3][2];
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 2; j ++)
      v[i][j] = V(j, vertices[i][0], vertices[i][1]);

  // check intersection
  double mu[3];
  bool succ = inverse_lerp_s2v2(v, mu);
  if (!succ) return false; // returns false if the cp is outside the triangle
  
  // lerp position
  double X[3][2];
  for (int i = 0; i < 3; i ++) 
    for (int j = 0; j < 2; j ++) 
      X[i][j] = vertices[i][j];
  lerp_s2v2(X, mu, cp.x);

  if (!scalar.empty()) { // lerp scalar value
    double values[3];
    for (int i = 0; i < 3; i ++)
      values[i] = scalar(vertices[i][0], vertices[i][1]);
    cp.scalar = lerp_s2(values, mu);
  }

  if (has_jacobian) {
    // derive jacobian
    double J[2][2]; // jacobian

    if (gradV.empty()) // jacobian is not given
      ftk::gradient_2dsimplex2_2(X, v, J);
    else { // lerp jacobian
      double Js[3][2][2];
      for (int i = 0; i < 3; i ++) 
        for (int j = 0; j < 2; j ++)
          for (int k = 0; k < 2; k ++)
            Js[i][j][k] = gradV(k, j, vertices[i][0], vertices[i][1]);
      lerp_s2m2x2(Js, mu, J);
    }
    cp.type = critical_point_type_2d(J, symmetric_jacobian);

    if (cp.type & type_filter) return true;
    else return false; // type mismatch
  } else return true;
}

#if FTK_HAVE_VTK
vtkSmartPointer<vtkPolyData> extract_critical_points_2d_regular_serial::get_results_vtk() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
  
  vtkIdType pid[1];
  for (const auto &cp : results) {
    double p[3] = {cp.x[0], cp.x[1], 0};
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  polyData->SetPoints(points);
  polyData->SetVerts(vertices);
 
  // point data for types
  vtkSmartPointer<vtkDoubleArray> types = vtkSmartPointer<vtkDoubleArray>::New();
  types->SetNumberOfValues(results.size());
  for (auto i = 0; i < results.size(); i ++) {
    types->SetValue(i, static_cast<double>(results[i].type));
  }
  types->SetName("type");
  polyData->GetPointData()->AddArray(types);
  
  // point data for scalars
  if (!scalar.empty()) {
    vtkSmartPointer<vtkDoubleArray> scalars = vtkSmartPointer<vtkDoubleArray>::New();
    scalars->SetNumberOfValues(results.size());
    for (auto i = 0; i < results.size(); i ++) {
      scalars->SetValue(i, static_cast<double>(results[i].scalar));
    }
    scalars->SetName("scalar");
    polyData->GetPointData()->AddArray(scalars);
  }

  return polyData;
}
#endif

}

#endif
