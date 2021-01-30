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
      float X[][4], float spin[][4]);
};

///////
inline void magnetic_vortex_tracker_3d_regular::simplex_values(
      const std::vector<std::vector<int>>& vertices,
      float X[][4], float spin[][3])
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    const auto &f = field_data_snapshots[iv];

    for (int k = 0; k < 3; k ++)
      spin[i][k] = f.spin(k, vertices[i][0], vertices[i][1], vertices[i][2]);
    spin[i][3] = 0;
    
    for (int j = 0; j < 3; j ++)
      X[i][j] = vertices[i][j];
    X[i][3] = vertices[i][3];
  }
}

inline bool magnetic_vortex_tracker_3d_regular::check_simplex(
    const element_t& e, feature_point_t& p)
{
  if (!e.valid(m)) return false; // check if the 2-simplex is valid
  const auto &vertices = e.vertices(m); // obtain the vertices of the simplex

  float X[3][4], // coordinates
        spin[3][4], 
        normal[4];
  simplex_values(vertices, X, spin);

  return false;

  return true;
}

inline void magnetic_vortex_tracker_3d_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
  
  auto get_relatetd_cels = [&](element_t e) {
    std::set<element_t> my_related_cells;
    
    auto tets = e.side_of(m);
    for (auto tet : tets) {
      if (tet.valid(m)) {
        auto pents = tet.side_of(m);
        for (auto pent : pents)
          if (pent.valid(m))
            my_related_cells.insert(pent);
      }
    }

    return my_related_cells;
  };

  auto func = [=](element_t e) {
    feature_point_t p;
    if (check_simplex(e, p)) {
      std::set<element_t> my_related_cells = get_relatetd_cels(e);

      {
        std::lock_guard<std::mutex> guard(mutex);
        intersections[e] = p;
        related_cells.insert(my_related_cells.begin(), my_related_cells.end());
      }
    }
  };
    
  element_for_ordinal(2, func);
  if (field_data_snapshots.size() >= 2) 
    element_for_interval(2, func);
}

} // namespace ftk

#endif
