#ifndef _FTK_PARTICLE_TRACER_HH
#define _FTK_PARTICLE_TRACER_HH

#include <ftk/config.hh>
#include <ftk/filters/tracker.hh>

namespace ftk {

struct particle_tracer : public virtual tracker
{
  particle_tracer(diy::mpi::communicator comm, int nd /*spatial dim*/) : tracker(comm), nd_(nd) {}
  virtual ~particle_tracer() {}

  virtual void initialize_particles(size_t n, const double *p, size_t stride=3, bool has_time=false);
  virtual void initialize_particles_at_grid_points() {};

  void update() {}
  void finalize() {}
  bool advance_timestep();

  void write_trajectories(const std::string& filename);

protected:
  // virtual bool eval_v(const double *x, double *v) = 0;

  int nd() const { return nd_; }
  double delta() const { return current_delta_t / nsteps_per_interval; }

protected:
  std::vector<feature_point_lite_t> particles;
  feature_curve_set_t trajectories;

  double current_t = 0.0, current_delta_t = 1.0;

  int nsteps_per_interval = 128;
  int nsteps_per_checkpoint = 16;

  const int nd_;
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

inline void particle_tracer::initialize_particles(size_t n, const double *buf, size_t stride, bool has_time)
{
  assert(stride >= nd());

  for (auto i = 0; i < n; i ++) {
    feature_curve_t curve;
    feature_point_t p;
   
    for (int k = 0; k < nd(); k ++)
      p.x[k] = buf[i*stride+k];

    if (has_time)
      p.t = buf[i*stride+nd()];

    curve.push_back(p);
    trajectories.add(curve);
  }
}

} // namespace ftk

#endif
