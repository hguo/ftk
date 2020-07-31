#ifndef _FTK_CRITICAL_POINT_TRACKER_2D_UNSTRUCTURED_HH
#define _FTK_CRITICAL_POINT_TRACKER_2D_UNSTRUCTURED_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/adjugate.hh>
#include <ftk/numeric/symmetric_matrix.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/geometry/curve2vtk.hh>
#include <ftk/filters/critical_point_tracker_regular.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/mesh/simplicial_unstructured_extruded_2d_mesh.hh>
#include <ftk/external/diy/serialization.hpp>

#if FTK_HAVE_GMP
#include <gmpxx.h>
#endif

#if FTK_HAVE_VTK
#include <vtkUnsignedIntArray.h>
#include <vtkVertex.h>
#endif

namespace ftk {

typedef critical_point_t<3, double> critical_point_2dt_t;

struct critical_point_tracker_2d_unstructured : public critical_point_tracker_regular
{
  // critical_point_tracker_2d_unstructured(const simplicial_unstructured_extruded_2d_mesh<>& m) : m(m) {}
  // critical_point_tracker_2d_unstructured() {}
  critical_point_tracker_2d_unstructured(const simplicial_unstructured_2d_mesh<>& m) : m(simplicial_unstructured_extruded_2d_mesh<>(m)) {}

  virtual ~critical_point_tracker_2d_unstructured() {};

  void initialize() {}
  void finalize();
  void reset() {}

  void update_timestep();

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_traced_critical_points_vtk() const;
  vtkSmartPointer<vtkPolyData> get_discrete_critical_points_vtk() const;
#endif

  void write_traced_critical_points_text(std::ostream& os) const {} // TODO
  void write_discrete_critical_points_text(std::ostream &os) const {} // TODO

  void push_scalar_field_snapshot(const ndarray<double>&) {} // TODO
  void push_vector_field_snapshot(const ndarray<double>&) {} // TODO

protected:
  bool check_simplex(int, critical_point_2dt_t& cp);

  template <typename T> void simplex_vectors(int n, int verts[], T v[][2]) const;
  void simplex_coordinates(int n, int verts[], double x[][3]) const;

protected:
  const simplicial_unstructured_extruded_2d_mesh<> m;
  
  std::map<int, critical_point_2dt_t> discrete_critical_points;
  std::vector<std::vector<critical_point_2dt_t>> traced_critical_points;
};

////////////////////////

template <typename T>
inline void critical_point_tracker_2d_unstructured::simplex_vectors(
    int n, int verts[], T v[][2]) const
{
  for (int i = 0; i < n; i ++) {
    const int iv = m.flat_vertex_time(verts[i]) == current_timestep ? 0 : 1;
    const int k = m.flat_vertex_id(verts[i]);
    for (int j = 0; j < 2; j ++) {
      v[i][j] = field_data_snapshots[iv].vector(j, k);
    }
  }
}

inline void critical_point_tracker_2d_unstructured::simplex_coordinates(
    int n, int verts[], double x[][3]) const
{
  for (int i = 0; i < n; i ++)
    m.get_coords(verts[i], x[i]);
}

inline bool critical_point_tracker_2d_unstructured::check_simplex(int i, critical_point_2dt_t& cp)
{
#if FTK_HAVE_GMP
  typedef mpf_class fp_t;
#else
  typedef ftk::fixed_point<> fp_t;
#endif
  
  int tri[3];
  m.get_simplex(2, i, tri); 

  double V[3][2], X[3][3];
  simplex_vectors<double>(3, tri, V);
  simplex_coordinates(3, tri, X);

  fp_t Vf[3][2];
  for (int k = 0; k < 3; k ++) 
    for (int j = 0; j < 2; j ++)
      Vf[k][j] = V[k][j];
   
  bool succ = ftk::robust_critical_point_in_simplex2(Vf, tri);
  if (!succ) return false;

  // ftk::print3x2("V", V);
  double mu[3], x[3];
  bool succ2 = ftk::inverse_lerp_s2v2(V, mu);
  // if (!succ2) return;
  ftk::lerp_s2v3(X, mu, cp.x);

#if 0
  double Js[3][2][2], H[2][2];
  for (int k = 0; k < 3; k ++) {
    int t = tri[k] >= m.n(0) ? 1 : 0;
    int v = tri[k] % m.n(0);
    for (int j = 0; j < 2; j ++)
      for (int i = 0; i < 2; i ++)
        Js[k][j][i] = J[t](i, j, v); 
  }
  ftk::lerp_s2m2x2(Js, mu, H);
  // ftk::print2x2("H", H);
  const int type = ftk::critical_point_type_2d(H, true);
#endif

  // fprintf(stderr, "mu=%f, %f, %f, x=%f, %f, %f, type=%d\n", 
  //     mu[0], mu[1], mu[2], x[0], x[1], x[2], type);

  cp.type = 0;
  return true;
  // cps.push_back(x[0]);
  // cps.push_back(x[1]);
  // cps.push_back(x[2]);
  // types.push_back(type);
}

inline void critical_point_tracker_2d_unstructured::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);

  auto func = [&](int i) {
    critical_point_2dt_t cp;
    if (check_simplex(i, cp)) {
      std::lock_guard<std::mutex> guard(mutex);
      discrete_critical_points[i] = cp;
    }
  };

  m.element_for_ordinal(2, current_timestep, func);
  if (field_data_snapshots.size() >= 2)
    m.element_for_interval(2, current_timestep, func);
}

inline void critical_point_tracker_2d_unstructured::finalize()
{
  fprintf(stderr, "finalizing...\n");
  // Convert connected components to geometries
  auto neighbors = [&](int f) {
    std::set<int> neighbors;
    const auto cells = m.side_of(2, f);
    fprintf(stderr, "face=%d\n", f);
    int vf[3];
    m.get_simplex(2, f, vf);
    fprintf(stderr, "face.simplex=%d, %d, %d\n", vf[0], vf[1], vf[2]);

    for (const auto c : cells) {
      fprintf(stderr, "--cell=%d\n", c);
      int vc[4];
      m.get_simplex(3, c, vc);
      fprintf(stderr, "--cell.simplex=%d, %d, %d, %d\n", vc[0], vc[1], vc[2], vc[3]);

      const auto elements = m.sides(3, c);
      for (const auto f1 : elements) {
        fprintf(stderr, "----face=%d\n", f1);
        neighbors.insert(f1);
      }
    }
    fprintf(stderr, "size_neighbors=%zu\n", neighbors.size());
    assert(neighbors.find(f) != neighbors.end());
    return neighbors;
  };

  std::set<int> elements;
  for (const auto &kv : discrete_critical_points)
    elements.insert(kv.first);
  auto connected_components = extract_connected_components<int, std::set<int>>(
      neighbors, elements);

  for (const auto &component : connected_components) {
    std::vector<std::vector<double>> mycurves;
    auto linear_graphs = ftk::connected_component_to_linear_components<int>(component, neighbors);
    fprintf(stderr, "size_component=%zu, size_linear_graph=%zu\n", component.size(), linear_graphs.size());
    for (int j = 0; j < linear_graphs.size(); j ++) {
      std::vector<critical_point_2dt_t> traj; 
      for (int k = 0; k < linear_graphs[j].size(); k ++)
        traj.push_back(discrete_critical_points[linear_graphs[j][k]]);
      traced_critical_points.emplace_back(traj);
    }
  }
  fprintf(stderr, "np=%zu, nc=%zu\n", discrete_critical_points.size(), connected_components.size());
}

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkPolyData> critical_point_tracker_2d_unstructured::get_discrete_critical_points_vtk() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
  
  vtkIdType pid[1];
  for (const auto &kv : discrete_critical_points) {
    const auto &cp = kv.second;
    double p[3] = {cp.x[0], cp.x[1], cp.x[2]};
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  polyData->SetPoints(points);
  polyData->SetVerts(vertices);
  return polyData;
}

inline vtkSmartPointer<vtkPolyData> critical_point_tracker_2d_unstructured::get_traced_critical_points_vtk() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> lines = vtkCellArray::New();
  vtkSmartPointer<vtkCellArray> verts = vtkCellArray::New();

  for (const auto &curve : traced_critical_points) {
    fprintf(stderr, "size=%zu\n", curve.size());
    for (auto i = 0; i < curve.size(); i ++) {
      double p[3] = {curve[i][0], curve[i][1], curve[i][2]};
      points->InsertNextPoint(p);
    }
  }

  size_t nv = 0;
  for (const auto &curve : traced_critical_points) {
    if (curve.size() < 2) { // isolated vertex
      vtkSmartPointer<vtkVertex> obj = vtkVertex::New();
      obj->GetPointIds()->SetNumberOfIds(curve.size());
      for (int i = 0; i < curve.size(); i ++)
        obj->GetPointIds()->SetId(i, i+nv);
      verts->InsertNextCell(obj);
    } else { // lines
      vtkSmartPointer<vtkPolyLine> obj = vtkPolyLine::New();
      obj->GetPointIds()->SetNumberOfIds(curve.size());
      for (int i = 0; i < curve.size(); i ++)
        obj->GetPointIds()->SetId(i, i+nv);
      lines->InsertNextCell(obj);
    }
    nv += curve.size();
  }
 
  polyData->SetPoints(points);
  polyData->SetLines(lines);
  polyData->SetVerts(verts);

  // point data for types
  if (1) { // if (type_filter) {
    vtkSmartPointer<vtkUnsignedIntArray> types = vtkSmartPointer<vtkUnsignedIntArray>::New();
    types->SetNumberOfValues(nv);
    size_t i = 0;
    for (const auto &curve : traced_critical_points) {
      for (auto j = 0; j < curve.size(); j ++)
        types->SetValue(i ++, curve[j].type);
    }
    types->SetName("type");
    polyData->GetPointData()->AddArray(types);
  }

  if (1) { // ids
    vtkSmartPointer<vtkUnsignedIntArray> ids = vtkSmartPointer<vtkUnsignedIntArray>::New();
    ids->SetNumberOfValues(nv);
    size_t i = 0;
    for (auto k = 0; k < traced_critical_points.size(); k ++)
      for (auto j = 0; j < traced_critical_points[k].size(); j ++)
        ids->SetValue(i ++, k);
    ids->SetName("id");
    polyData->GetPointData()->AddArray(ids);

  }

  // point data for scalars
  // if (has_scalar_field) {
  if (1) { // scalar is 0 if no scalar field available
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
