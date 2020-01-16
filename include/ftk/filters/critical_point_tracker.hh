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

#if FTK_HAVE_VTK
  virtual vtkSmartPointer<vtkPolyData> get_traced_critical_points_vtk() const = 0;
  virtual vtkSmartPointer<vtkPolyData> get_discrete_critical_points_vtk() const = 0;
#endif
  void write_traced_critical_points_vtk(const std::string& filename);
  void write_discrete_critical_points_vtk(const std::string& filename);
};

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

}

#endif
