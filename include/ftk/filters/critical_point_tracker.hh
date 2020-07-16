#ifndef _FTK_CRITICAL_POINT_TRACKER_HH
#define _FTK_CRITICAL_POINT_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/filters/filter.hh>
#include <ftk/filters/critical_point.hh>
#include <ftk/geometry/points2vtk.hh>

namespace ftk {

struct critical_point_tracker : public filter {
  critical_point_tracker(int argc, char **argv) : filter(argc, argv) {}
  critical_point_tracker() {}

  virtual void update() {}; // TODO
  void reset() {field_data_snapshots.clear();}

#if FTK_HAVE_VTK
  virtual vtkSmartPointer<vtkPolyData> get_traced_critical_points_vtk() const = 0;
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

protected:
  std::deque<field_data_snapshot_t> field_data_snapshots;
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
  if (comm.rank() == 0) {
    auto poly = get_traced_critical_points_vtk();
    write_vtp(filename, poly);
  }
}

inline void critical_point_tracker::write_discrete_critical_points_vtk(const std::string& filename)
{
  if (comm.rank() == 0) {
    auto poly = get_discrete_critical_points_vtk();
    write_vtp(filename, poly);
  }
}
#else
inline void critical_point_tracker::write_traced_critical_points_vtk(const std::string& filename)
{
  if (comm.rank() == 0)
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}

inline void critical_point_tracker::write_discrete_critical_points_vtk(const std::string& filename)
{
  if (comm.rank() == 0)
    fprintf(stderr, "[FTK] fatal: FTK not compiled with VTK.\n");
}
#endif

inline void critical_point_tracker::write_traced_critical_points_text(const std::string& filename)
{
  if (comm.rank() == 0) {
    std::ofstream out(filename);
    write_traced_critical_points_text(out);
    out.close();
  }
}

inline void critical_point_tracker::write_discrete_critical_points_text(const std::string& filename)
{
  if (comm.rank() == 0) {
    std::ofstream out(filename);
    write_discrete_critical_points_text(out);
    out.close();
  }
}

}

#endif
