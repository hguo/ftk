#ifndef _FTK_MAGNETIC_VORTEX_TRACKER_3D_REGULAR_HH
#define _FTK_MAGNETIC_VORTEX_TRACKER_3D_REGULAR_HH

#include <ftk/config.hh>
#include <ftk/filters/magnetic_vortex_tracker.hh>
#include <ftk/filters/regular_tracker.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/write_polydata.hh>
#include <ftk/ndarray/writer.hh>

namespace ftk {

struct magnetic_vortex_tracker_3d_regular : public virtual magnetic_vortex_tracker, public virtual regular_tracker
{
  magnetic_vortex_tracker_3d_regular(diy::mpi::communicator comm) : tdgl_vortex_tracker(comm), regular_tracker(comm, 3), tracker(comm) {}
  virtual ~magnetic_vortex_tracker_3d_regular();

  void finalize();
  void reset();

  void update_timestep();

protected:
  typedef simplicial_regular_mesh_element element_t;
  
  std::map<element_t, feature_point_t> intersections;
  std::set<element_t> related_cells;

  feature_surface_t surfaces;

protected:
  bool check_simplex(const element_t& s, feature_point_t& cp);
  
  void simplex_values(
      const std::vector<std::vector<int>>& vertices,
      float spin[][3]);
};

}

#endif
