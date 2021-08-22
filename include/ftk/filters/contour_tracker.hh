#ifndef _FTK_CONTOUR_TRACKER_HH
#define _FTK_CONTOUR_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/algorithms/cca.hh>
// #include <ftk/filters/contour.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/features/feature_point.hh>
#include <ftk/features/feature_surface.hh>
#include <ftk/features/feature_volume.hh>
#include <ftk/geometry/write_polydata.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/utils/gather.hh>
#include <iomanip>

namespace ftk {

struct contour_tracker : public virtual tracker {
  contour_tracker(diy::mpi::communicator comm) : tracker(comm) {}

  virtual void update() {}; 
  void reset() {
    field_data_snapshots.clear();
    // traced_contours.clear();
  }
  
  void set_scalar_components(const std::vector<std::string>& c);
  int get_num_scalar_components() const {return scalar_components.size();}

  double get_threshold() const { return threshold; }
  void set_threshold(double t) {threshold = t;}

public:
  virtual bool advance_timestep();

public: // inputs
  bool pop_field_data_snapshot();
  virtual void push_field_data_snapshot(const ndarray<double> &scalar);

public:
  virtual std::vector<feature_point_t> get_intersections() const = 0;

  void write_intersections(const std::string& filenames) const;
  void write_intersections_vtp(const std::string& filenames) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_intersections_vtp() const;
#endif

  virtual void write_isovolume_vtu(const std::string& filename) const = 0;
  virtual void write_sliced_vtu(const std::string& pattern) const {}
  virtual void write_sliced_vtp(const std::string& pattern) const {}

protected:
  // virtual int cpdims() const = 0;

protected:
  struct field_data_snapshot_t {
    ndarray<double> scalar, gradient;
  };
  std::deque<field_data_snapshot_t> field_data_snapshots;
  
  double threshold = 0.0;
  
  // scalar components
  std::vector<std::string> scalar_components = {"scalar"};
};

///////

inline void contour_tracker::push_field_data_snapshot(const ndarray<double>& scalar)
{
  field_data_snapshot_t snapshot;
  snapshot.scalar = scalar;
  if (scalar.nd() == 2) 
    snapshot.gradient = gradient2D(scalar);
  else 
    snapshot.gradient = gradient3D(scalar);

  field_data_snapshots.emplace_back(snapshot);
}

inline bool contour_tracker::pop_field_data_snapshot()
{
  if (field_data_snapshots.size() > 0) {
    field_data_snapshots.pop_front();
    return true;
  } else return false;
}

inline bool contour_tracker::advance_timestep()
{
  update_timestep();
  pop_field_data_snapshot();

  current_timestep ++;
  return field_data_snapshots.size() > 0;
}

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkPolyData> contour_tracker::get_intersections_vtp() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
 
  vtkIdType pid[1];
  for (const auto &cp : get_intersections()) {
    double p[3] = {cp.x[0], cp.x[1], cp.x[2]}; 
    // if (cpdims() == 2) p[2] = cp.t;
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  polyData->SetPoints(points);
  polyData->SetVerts(vertices);
  
  vtkSmartPointer<vtkDoubleArray> time_array = vtkSmartPointer<vtkDoubleArray>::New();
  time_array->SetNumberOfValues(get_intersections().size());
  int i = 0;
  for (const auto &cp : get_intersections())
    time_array->SetValue(i++, cp.t);
  time_array->SetName("time");
  polyData->GetPointData()->AddArray(time_array);

  return polyData;
}

inline void contour_tracker::write_intersections_vtp(const std::string& filename) const
{
  if (comm.rank() == get_root_proc()) {
    auto poly = get_intersections_vtp();
    write_polydata(filename, poly);
  }
}
#else
inline void contour_tracker::write_intersections_vtp(const std::string& filename) const
{
  if (is_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}
#endif

}

#endif
