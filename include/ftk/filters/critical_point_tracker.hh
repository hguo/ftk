#ifndef _FTK_CRITICAL_POINT_TRACKER_HH
#define _FTK_CRITICAL_POINT_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/filters/filter.hh>
#include <ftk/filters/critical_point.hh>
#include <ftk/geometry/points2vtk.hh>

#if FTK_HAVE_VTK
#include <vtkUnsignedIntArray.h>
#include <vtkVertex.h>
#endif

namespace ftk {

struct critical_point_tracker : public filter {
  critical_point_tracker() {}

  virtual void update() {}; // TODO
  void reset() {
    field_data_snapshots.clear();
    traced_critical_points.clear();
  }

#if FTK_HAVE_VTK
  virtual vtkSmartPointer<vtkPolyData> get_traced_critical_points_vtk() const;
  virtual vtkSmartPointer<vtkPolyData> get_discrete_critical_points_vtk() const = 0;
#endif
  void write_traced_critical_points_vtk(const std::string& filename);
  void write_discrete_critical_points_vtk(const std::string& filename);

  virtual void write_traced_critical_points_text(std::ostream& os) const = 0;
  virtual void write_discrete_critical_points_text(std::ostream &os) const = 0;

  void write_traced_critical_points_text(const std::string& filename);
  void write_discrete_critical_points_text(const std::string& filename);

  struct field_data_snapshot_t {
    ndarray<double> scalar, vector, jacobian;
  };

  bool pop_field_data_snapshot();
  virtual void push_field_data_snapshot(
      const ndarray<double> &scalar, 
      const ndarray<double> &vector,
      const ndarray<double> &jacobian);
  virtual void push_scalar_field_snapshot(const ndarray<double> &scalar); // push scalar only

  virtual void push_field_data_spacetime(
      const ndarray<double> &scalars, 
      const ndarray<double> &vectors,
      const ndarray<double> &jacobians);
  void push_scalar_field_spacetime(const ndarray<double>& scalars);

  const std::vector<std::vector<critical_point_t>>& get_traced_critical_points() {return traced_critical_points;}

protected:
  std::deque<field_data_snapshot_t> field_data_snapshots;
  
  std::vector<std::vector<critical_point_t>> traced_critical_points;
};

///////

inline void critical_point_tracker::push_field_data_snapshot(
    const ndarray<double>& scalar,
    const ndarray<double>& vector,
    const ndarray<double>& jacobian)
{
  field_data_snapshot_t snapshot;
  snapshot.scalar = scalar;
  snapshot.vector = vector;
  snapshot.jacobian = jacobian;

  field_data_snapshots.emplace_back(snapshot);
}

inline void critical_point_tracker::push_scalar_field_snapshot(const ndarray<double>& scalar)
{
  field_data_snapshot_t snapshot;
  snapshot.scalar = scalar;

  field_data_snapshots.emplace_back(snapshot);
}

inline void critical_point_tracker::push_field_data_spacetime(
    const ndarray<double>& scalars,
    const ndarray<double>& vectors,
    const ndarray<double>& jacobians)
{
  for (size_t t = 0; t < scalars.shape(scalars.nd()-1); t ++) {
    auto scalar = scalars.slice_time(t);
    auto vector = vectors.slice_time(t);
    auto jacobian = jacobians.slice_time(t);

    push_field_data_snapshot(scalar, vector, jacobian);
  }
}

inline void critical_point_tracker::push_scalar_field_spacetime(const ndarray<double>& scalars)
{
  for (size_t t = 0; t < scalars.shape(scalars.nd()-1); t ++)
    push_scalar_field_snapshot( scalars.slice_time(t) );
}


inline bool critical_point_tracker::pop_field_data_snapshot()
{
  if (field_data_snapshots.size() > 0) {
    field_data_snapshots.pop_front();
    return true;
  } else return false;
}

//////
#if FTK_HAVE_VTK
inline void critical_point_tracker::write_traced_critical_points_vtk(const std::string& filename)
{
  if (comm.rank() == get_root_proc()) {
    auto poly = get_traced_critical_points_vtk();
    write_vtp(filename, poly);
  }
}

inline void critical_point_tracker::write_discrete_critical_points_vtk(const std::string& filename)
{
  if (comm.rank() == get_root_proc()) {
    auto poly = get_discrete_critical_points_vtk();
    write_vtp(filename, poly);
  }
}
#else
inline void critical_point_tracker::write_traced_critical_points_vtk(const std::string& filename)
{
  if (comm.rank() == get_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}

inline void critical_point_tracker::write_discrete_critical_points_vtk(const std::string& filename)
{
  if (comm.rank() == get_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}
#endif

inline void critical_point_tracker::write_traced_critical_points_text(const std::string& filename)
{
  if (comm.rank() == get_root_proc()) {
    std::ofstream out(filename);
    write_traced_critical_points_text(out);
    out.close();
  }
}

inline void critical_point_tracker::write_discrete_critical_points_text(const std::string& filename)
{
  if (comm.rank() == get_root_proc()) {
    std::ofstream out(filename);
    write_discrete_critical_points_text(out);
    out.close();
  }
}

inline vtkSmartPointer<vtkPolyData> critical_point_tracker::get_traced_critical_points_vtk() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> lines = vtkCellArray::New();
  vtkSmartPointer<vtkCellArray> verts = vtkCellArray::New();

  for (const auto &curve : traced_critical_points)
    for (auto i = 0; i < curve.size(); i ++) {
      double p[3] = {curve[i][0], curve[i][1], curve[i][2]};
      points->InsertNextPoint(p);
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
        scalars->SetValue(i ++, curve[j].scalar[0]);
    }
    scalars->SetName("scalar");
    polyData->GetPointData()->AddArray(scalars);
  }

  return polyData;
}

}

#endif
