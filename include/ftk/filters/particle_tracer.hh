#ifndef _FTK_PARTICLE_TRACER_HH
#define _FTK_PARTICLE_TRACER_HH

#include <ftk/config.hh>
#include <ftk/filters/tracker.hh>

namespace ftk {

struct particle_tracer : public virtual tracker
{
  particle_tracer(diy::mpi::communicator comm) : tracker(comm) {}
  virtual ~particle_tracer() {}

  virtual void set_particles(const ndarray<double> &p) {} //  particles = p; }
  virtual void initialize_particles_at_grid_points() {};

  void update() {}
  void finalize() {}
  bool advance_timestep();

  void write_trajectories(const std::string& filename);

protected:
  // virtual bool eval_v(const double *x, double *v) = 0;

protected:
  std::vector<feature_point_lite_t> particles;
  feature_curve_set_t trajectories;

  double current_t = 0.0, current_delta_t = 1.0;
};

//// 
inline bool particle_tracer::advance_timestep()
{
  update_timestep();
  pop_field_data_snapshot();

  current_timestep ++;
  return snapshots.size() > 0;
}

inline void particle_tracer::write_trajectories(const std::string& filename)
{
  // if (comm.rank() == 0)
#if FTK_HAVE_VTK
  auto poly = trajectories.to_vtp();
  write_polydata(filename, poly);
#endif
}

}

#endif
