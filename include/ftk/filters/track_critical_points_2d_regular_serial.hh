#ifndef _FTK_TRACK_CRITICAL_POINTS_2D_REGULAR_SERIAL_HH
#define _FTK_TRACK_CRITICAL_POINTS_2D_REGULAR_SERIAL_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/geometry/curve2vtk.hh>
#include <ftk/filters/extract_critical_points_2d_regular_serial.hh>
#include <ftk/ndarray.hh>
#include <ftk/hypermesh/regular_simplex_mesh.hh>

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
  void set_input_scalar_field(const ndarray<double>& scalar_) {scalar = scalar_;}

  void set_input_vector_field(const double *p, size_t W, size_t H, size_t T);
  void set_input_vector_field(const ndarray<double>& V_) {V = V_;}

  void set_input_jacobian_field(const double *p, size_t W, size_t H, size_t T); 
  void set_input_jacobian_field(const ndarray<double> &J) {gradV = J;}
  
  void set_lb_ub(const std::vector<int>& lb, const std::vector<int>& ub) {m.set_lb_ub(lb, ub);}
  void set_type_filter(unsigned int mask = 0xffffffff) {type_filter = mask;}

  void get_results();

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_results_vtk() const;
#endif

protected:
  ndarray<double> scalar, V, gradV;
  regular_simplex_mesh m;
  
  unsigned int type_filter = 0xffffffff;
  bool has_jacobian = false;
  bool symmetric_jacobian = false;

  typedef regular_simplex_mesh_element element_t;
  
  std::map<element_t, critical_point_2dt_t> discrete_critical_points;
  std::vector<std::set<element_t>> connected_components;
  std::vector<std::vector<critical_point_2dt_t>> traced_critical_points;

protected:
  bool check_simplex(const element_t& s, critical_point_2dt_t& cp);
  void trace_connected_components();
};


////////////////////
void track_critical_points_2d_regular_serial::execute()
{
  // initializing vector fields
  if (!scalar.empty()) {
    if (V.empty()) V = gradient2Dt(scalar);
    if (gradV.empty()) gradV = jacobian2Dt(V);
    has_jacobian = true;
    symmetric_jacobian = true;
  }

  // initializing bounds
  if (m.lb() == m.ub()) {
    if (!scalar.empty())
      m.set_lb_ub({2, 2, 0}, {static_cast<int>(V.dim(1)-3), static_cast<int>(V.dim(2)-3), static_cast<int>(V.dim(3)-1)});
    else
      m.set_lb_ub({0, 0, 0}, {static_cast<int>(V.dim(1)-1), static_cast<int>(V.dim(2)-1), static_cast<int>(V.dim(3)-1)});
  }

  // scan 2-simplices
  fprintf(stderr, "tracking 2D critical points...\n");
  m.element_for(2, [=](regular_simplex_mesh_element e) {
      critical_point_2dt_t cp;
      if (check_simplex(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        discrete_critical_points[e] = cp;
      }
    }); 

  // scan 3-simplices to get connected components
  fprintf(stderr, "trace intersections...\n");
  union_find<element_t> uf;
  for (const auto &kv : discrete_critical_points) 
    uf.add(kv.first);

  m.element_for(3, [&](const regular_simplex_mesh_element& f) {
    const auto sides = f.sides();
    std::set<element_t> intersected_sides;

    for (const auto& side : sides)
      if (uf.has(side)) 
        intersected_sides.insert(side);

    if (intersected_sides.size() > 1) {// the size of intersected_size should be 0 or 2
      for (auto it = std::next(intersected_sides.begin(), 1); it != intersected_sides.end(); it ++) {
        std::lock_guard<std::mutex> guard(mutex);
        uf.unite(*intersected_sides.begin(), *it);
      }
    }
  });
  uf.get_sets(connected_components);

  // convert connected components to traced critical points
  fprintf(stderr, "tracing critical points...\n");
  trace_connected_components();
}

void track_critical_points_2d_regular_serial::trace_connected_components()
{
  // Convert connected components to geometries
  auto neighbors = [](element_t f) {
    std::set<element_t> neighbors;
    const auto cells = f.side_of();
    for (const auto c : cells) {
      const auto elements = c.sides();
      for (const auto f1 : elements)
        neighbors.insert(f1);
    }
    return neighbors;
  };

  for (const auto &component : connected_components) {
    std::vector<std::vector<double>> mycurves;
    auto linear_graphs = ftk::connected_component_to_linear_components<element_t>(component, neighbors);
    for (int j = 0; j < linear_graphs.size(); j ++) {
      std::vector<critical_point_2dt_t> traj; 
      for (int k = 0; k < linear_graphs[j].size(); k ++)
        traj.push_back(discrete_critical_points[linear_graphs[j][k]]);
      traced_critical_points.emplace_back(traj);
    }
  }
}

bool track_critical_points_2d_regular_serial::check_simplex(
    const regular_simplex_mesh_element& e,
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

  if (!scalar.empty()) {
    double values[3];
    for (int i = 0; i < 3; i ++)
      values[i] = scalar(vertices[i][0], vertices[i][1], vertices[i][2]);
    cp.scalar = lerp_s2(values, mu);
  }

  if (has_jacobian) {
    // derive jacobian
    double J[2][2] = {0};
    if (gradV.empty()) {
      // TODO: jacobian is not given
    } else { // lerp jacobian
      double Js[3][2][2];
      for (int i = 0; i < 3; i ++)
        for (int j = 0; j < 2; j ++)
          for (int k = 0; k < 2; k ++) 
            Js[i][j][k] = gradV(k, j, vertices[i][0], vertices[i][1], vertices[i][2]);
      lerp_s2m2x2(Js, mu, J);
    }
    cp.type = critical_point_type_2d(J, symmetric_jacobian);
    if (cp.type & type_filter) return true;
    else return false;
  } else return true;
}

#if FTK_HAVE_VTK
vtkSmartPointer<vtkPolyData> track_critical_points_2d_regular_serial::get_results_vtk() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> cells = vtkCellArray::New();

  for (const auto &curve : traced_critical_points)
    for (auto i = 0; i < curve.size(); i ++) {
      double p[3] = {curve[i][0], curve[i][1], curve[i][2]};
      points->InsertNextPoint(p);
    }

  size_t nv = 0;
  for (const auto &curve : traced_critical_points) {
    vtkSmartPointer<vtkPolyLine> polyLine = vtkPolyLine::New();
    polyLine->GetPointIds()->SetNumberOfIds(curve.size());
    for (int i = 0; i < curve.size(); i ++)
      polyLine->GetPointIds()->SetId(i, i+nv);

    cells->InsertNextCell(polyLine);
    nv += curve.size();
  }
  
  polyData->SetPoints(points);
  polyData->SetLines(cells);

  // point data for types
  if (type_filter) {
    vtkSmartPointer<vtkDoubleArray> types = vtkSmartPointer<vtkDoubleArray>::New();
    types->SetNumberOfValues(nv);
    size_t i = 0;
    for (const auto &curve : traced_critical_points) {
      for (auto j = 0; j < curve.size(); j ++)
        types->SetValue(i ++, curve[j].type);
    }
    types->SetName("type");
    polyData->GetPointData()->AddArray(types);
  }

  // point data for scalars
  if (!scalar.empty()) {
    vtkSmartPointer<vtkDoubleArray> scalars = vtkSmartPointer<vtkDoubleArray>::New();
    scalars->SetNumberOfValues(nv);
    size_t i = 0;
    for (const auto &curve : traced_critical_points) {
      for (auto j = 0; j < curve.size(); j ++)
        scalars->SetValue(i ++, curve[j].scalar);
    }
    scalars->SetName("scalar");
    polyData->GetPointData()->AddArray(scalars);
  }

  return polyData;
}
#endif

}

#endif
