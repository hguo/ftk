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

protected:
  virtual bool eval_v(const double *x, double *v) = 0;

protected:
  std::vector<feature_point_lite_t> particles;

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

}

#endif
