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

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_results_vtk() const;
#endif

protected:
  hypermesh::ndarray<double> scalar, V, gradV;
  hypermesh::regular_simplex_mesh m;
  
  unsigned int type_filter = 0xffffffff;
  unsigned int jacobian_mode = JACOBIAN_NONE;

  std::map<hypermesh::regular_simplex_mesh_element, critical_point_2dt_t> intersections;
  std::vector<std::vector<double>> trajectories; // the output trajectories

protected:
  bool check_simplex(const hypermesh::regular_simplex_mesh_element& s, critical_point_2dt_t& cp);
  void trace_intersections();
  void extract_connected_components(std::vector<std::set<hypermesh::regular_simplex_mesh_element>>& components);
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

  fprintf(stderr, "trace intersections...\n");
  trace_intersections();
}

void track_critical_points_2d_regular_serial::extract_connected_components(std::vector<std::set<hypermesh::regular_simplex_mesh_element>>& components)
{
  typedef hypermesh::regular_simplex_mesh_element element_t;

  // Initialization
  ftk::union_find<std::string> uf; 
  std::map<std::string, element_t> id2ele; 
  for (const auto &f : intersections) {
    std::string eid = f.first.to_string(); 
    uf.add(eid); 
    id2ele.insert(std::make_pair(eid, f.first));
  }

  // Connected Component Labeling by using union-find. 
  m.element_for(3, [&](const hypermesh::regular_simplex_mesh_element& f) {
    const auto elements = f.sides();
    std::set<std::string> features; 

    for (const auto& ele : elements) {
      std::string eid = ele.to_string(); 

      if(uf.has(eid)) {
        features.insert(eid); 
      }
    }

    if(features.size()  > 1) {
      for(std::set<std::string>::iterator ite_i = std::next(features.begin(), 1); ite_i != features.end(); ++ite_i) {
        std::lock_guard<std::mutex> guard(mutex); // Use a lock for thread-save. 
        uf.unite(*(features.begin()), *ite_i); 
      }
    }
  }); 


  // the output sets of connected elements
  std::vector<std::set<std::string>> connected_components_str; // connected components 
  
  // Get disjoint sets of element IDs
  uf.get_sets(connected_components_str);

  // Convert element IDs to elements
  for(auto& comp_str : connected_components_str) {
    std::set<element_t>& comp = components.emplace_back(); 
    for(auto& ele_id : comp_str) {
      comp.insert(id2ele.find(ele_id)->second); 
    }
  }
}

void track_critical_points_2d_regular_serial::trace_intersections()
{
  typedef hypermesh::regular_simplex_mesh_element element_t; 

  std::vector<std::set<element_t>> cc; // connected components 
  extract_connected_components(cc);

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

  fprintf(stderr, "generating curves...\n");
#if 1
  for (int i = 0; i < cc.size(); i ++) {
    std::vector<std::vector<double>> mycurves;
    auto linear_graphs = ftk::connected_component_to_linear_components<element_t>(cc[i], neighbors);
    for (int j = 0; j < linear_graphs.size(); j ++) {
      std::vector<double> mycurve, mycolors;
      // float max_value = std::numeric_limits<float>::min();
      for (int k = 0; k < linear_graphs[j].size(); k ++) {
        auto p = intersections[linear_graphs[j][k]];
        mycurve.push_back(p.x[0]); //  / (DW-1));
        mycurve.push_back(p.x[1]); //  / (DH-1));
        mycurve.push_back(p.x[2]); //  / (DT-1));
        mycurve.push_back(p.scalar);
        // max_value = std::max(max_value, p.scalar);
      }
      trajectories.emplace_back(mycurve);
    }
  }
#endif
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
}

#if FTK_HAVE_VTK
vtkSmartPointer<vtkPolyData> track_critical_points_2d_regular_serial::get_results_vtk() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
  
  vtkIdType pid[1];
  for (const auto &kv : intersections) {
    const auto &cp = kv.second;
    double p[3] = {cp.x[0], cp.x[1], cp.x[2]};
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  polyData->SetPoints(points);
  polyData->SetVerts(vertices);
 
  return polyData;
}
#endif

}

#endif
