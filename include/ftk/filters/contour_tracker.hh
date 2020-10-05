#ifndef _FTK_CRITICAL_POINT_TRACKER_HH
#define _FTK_CRITICAL_POINT_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/algorithms/cca.hh>
// #include <ftk/filters/contour.hh>
#include <ftk/filters/filter.hh>
#include <ftk/filters/feature_point.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/external/diy/serialization.hpp>
#include <ftk/external/diy-ext/gather.hh>
#include <iomanip>

namespace ftk {

struct contour_tracker : public filter {
  contour_tracker() {}

  virtual void update() {}; 
  void reset() {
    field_data_snapshots.clear();
    // traced_contours.clear();
  }
  
  void set_scalar_components(const std::vector<std::string>& c);
  int get_num_scalar_components() const {return scalar_components.size();}

  void set_threshold(double t) {threshold = t;}

public:
  virtual void initialize() = 0;
  virtual void finalize() = 0;

  virtual bool advance_timestep();
  virtual void update_timestep() = 0;

public: // inputs
  bool pop_field_data_snapshot();
  virtual void push_field_data_snapshot(const ndarray<double> &scalar);

  virtual void set_current_timestep(int t) {current_timestep = t;}
  int get_current_timestep() const {return current_timestep;}

  void set_end_timestep(int t) {end_timestep = t;}

public:
  virtual std::vector<feature_point_t> get_intersections() const = 0;

  void write_intersections(const std::string& filenames) const;
  void write_intersections_vtk(const std::string& filenames) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkPolyData> get_intersections_vtk() const;
#endif

  virtual void write_trajectories_vtk(const std::string& filename) const = 0;

protected:
  virtual int cpdims() const = 0;

protected:
  struct field_data_snapshot_t {
    ndarray<double> scalar;
  };
  std::deque<field_data_snapshot_t> field_data_snapshots;
  
  int current_timestep = 0;
  int start_timestep = 0, 
      end_timestep = std::numeric_limits<int>::max();

  double threshold = 0.0;
  
  // scalar components
  std::vector<std::string> scalar_components = {"scalar"};
};

///////

inline void contour_tracker::push_field_data_snapshot(const ndarray<double>& scalar)
{
  field_data_snapshot_t snapshot;
  snapshot.scalar = scalar;

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
inline vtkSmartPointer<vtkPolyData> contour_tracker::get_intersections_vtk() const
{
  vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
  
  vtkIdType pid[1];
  
  for (const auto &cp : get_intersections()) {
    // double p[3] = {cp.x[0], cp.x[1], cp.x[2]}; // TODO: time
    double p[3] = {cp.x[0], cp.x[1], cp.t};
    if (cpdims() == 3) p[2] = cp.x[2];
    pid[0] = points->InsertNextPoint(p);
    vertices->InsertNextCell(1, pid);
  }

  polyData->SetPoints(points);
  polyData->SetVerts(vertices);

  return polyData;
}

inline void contour_tracker::write_intersections_vtk(const std::string& filename) const
{
  if (comm.rank() == get_root_proc()) {
    auto poly = get_intersections_vtk();
    write_vtp(filename, poly);
  }
}
#else
inline void contour_tracker::write_intersections_vtk(const std::string& filename) const
{
  if (is_root_proc())
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}
#endif

}

#endif
