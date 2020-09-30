#ifndef _FTK_CONTOUR_TRACKER_2D_REGULAR_HH
#define _FTK_CONTOUR_TRACKER_2D_REGULAR_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/adjugate.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/geometry/curve2vtk.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/filters/contour_tracker_regular.hh>
#include <ftk/external/diy/serialization.hpp>
#include <ftk/external/diy-ext/gather.hh>

#if FTK_HAVE_GMP
#include <gmpxx.h>
#endif

namespace ftk {

struct contour_tracker_2d_regular : public contour_tracker_regular {
  contour_tracker_2d_regular() : contour_tracker_regular(2) {}
  virtual ~contour_tracker_2d_regular() {}

  int cpdims() const { return 2; }

  void initialize();
  void finalize();
  void reset();

  void update_timestep();

protected:
  typedef simplicial_regular_mesh_element element_t;
  
protected:
  bool check_simplex(const element_t& s, feature_point_t& cp);
  void trace_intersections();
  void trace_connected_components();

  virtual void simplex_coordinates(const std::vector<std::vector<int>>& vertices, double X[][3]) const;
  virtual void simplex_scalars(const std::vector<std::vector<int>>& vertices, double values[]) const;
};


////////////////////
inline void contour_tracker_2d_regular::initialize()
{
  // initializing bounds
  m.set_lb_ub({
      static_cast<int>(domain.start(0)),
      static_cast<int>(domain.start(1)),
      start_timestep
    }, {
      static_cast<int>(domain.size(0)),
      static_cast<int>(domain.size(1)),
      end_timestep
    });
}

inline void contour_tracker_2d_regular::finalize()
{
  diy::mpi::gather(comm, intersections, intersections, get_root_proc());

  if (comm.rank() == get_root_proc()) {
    // fprintf(stderr, "finalizing...\n");
    // traced_critical_points.add( trace_critical_points_offline<element_t>(discrete_critical_points, 
    //     [&](element_t f) {
    //       std::set<element_t> neighbors;
    //       const auto cells = f.side_of(m);
    //       for (const auto c : cells) {
    //         const auto elements = c.sides(m);
    //         for (const auto f1 : elements)
    //           neighbors.insert(f1);
    //       }
    //       return neighbors;
    // }));
  }
}

inline void contour_tracker_2d_regular::reset()
{
  current_timestep = 0;

  field_data_snapshots.clear();
  intersections.clear();

  contour_tracker::reset();
}

inline void contour_tracker_2d_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
}

inline void contour_tracker_2d_regular::simplex_coordinates(
    const std::vector<std::vector<int>>& vertices, double X[][3]) const
{
  for (int i = 0; i < vertices.size(); i ++)
    for (int j = 0; j < 3; j ++)
      X[i][j] = vertices[i][j];
}

inline void contour_tracker_2d_regular::simplex_scalars(
    const std::vector<std::vector<int>>& vertices, double values[]) const
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][2] == current_timestep ? 0 : 1;
    values[i] = field_data_snapshots[iv].scalar(
        vertices[i][0],
        vertices[i][1]);
  }
}

}

#endif
