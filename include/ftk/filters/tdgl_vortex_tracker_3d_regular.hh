#ifndef _FTK_TDGL_TRACKER_3D_REGULAR_HH
#define _FTK_TDGL_TRACKER_3D_REGULAR_HH

#include <ftk/ftk_config.hh>
#include <ftk/filters/tdgl_vortex_tracker.hh>

namespace ftk {

struct tdgl_vortex_tracker_3d_regular : public tdgl_vortex_tracker
{
  tdgl_vortex_tracker_3d_regular() : m(4) {}
  virtual ~tdgl_vortex_tracker_3d_regular() {}
  
  void set_domain(const lattice& l) {domain = l;} // spatial domain
  void set_array_domain(const lattice& l) {array_domain = l;}

  void initialize();
  void finalize();
  void reset();

  void update_timestep();
 
protected:
  void simplex_values(
      const std::vector<std::vector<int>>& vertices,
      float X[][4],
      float rho[], float phi[],
      float re[], float im[]);

protected:
  simplicial_regular_mesh m;
  typedef simplicial_regular_mesh_element element_t;
  
  std::map<element_t, feature_point_t> intersections;
  std::set<element_t> related_cells;

protected: // config
  lattice domain, array_domain;
};


////////////////////
inline void tdgl_vortex_tracker_3d_regular::initialize()
{
  // initializing bounds
  m.set_lb_ub({
      static_cast<int>(domain.start(0)),
      static_cast<int>(domain.start(1)),
      static_cast<int>(domain.start(2)),
      start_timestep
    }, {
      static_cast<int>(domain.size(0)),
      static_cast<int>(domain.size(1)),
      static_cast<int>(domain.size(2)),
      end_timestep
    });
}

inline void tdgl_vortex_tracker_3d_regular::finalize()
{
  // TODO
}

inline void tdgl_vortex_tracker_3d_regular::reset()
{
  current_timestep = 0;

  field_data_snapshots.clear();
  intersections.clear();

  tdgl_vortex_tracker::reset();
}

inline void tdgl_vortex_tracker_3d_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
}
  
inline void tdgl_vortex_tracker_3d_regular::simplex_values(
      const std::vector<std::vector<int>>& vertices,
      float X[][4],
      float rho[], float phi[],
      float re[], float im[])
{
  for (int i = 0; i < vertices.size(); i ++) {
    const int iv = vertices[i][3] == current_timestep ? 0 : 1;
    const auto &f = field_data_snapshots[iv];

    rho[i] = f.rho(vertices[i][0], vertices[i][1], vertices[i][2]);
    phi[i] = f.phi(vertices[i][0], vertices[i][1], vertices[i][2]);
    re[i] = f.re(vertices[i][0], vertices[i][1], vertices[i][2]);
    im[i] = f.im(vertices[i][0], vertices[i][1], vertices[i][2]);
    
    for (int j = 0; j < 4; j ++)
      X[i][j] = vertices[i][j];
  }
}

}

#endif
