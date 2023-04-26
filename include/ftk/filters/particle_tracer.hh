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
  void update_timestep();
  bool advance_timestep();

  void write_trajectories(const std::string& filename);

protected:
  // virtual bool eval_v(const double *x, double *v) = 0;
  virtual bool eval_v(std::shared_ptr<ndarray<double>> V0,
      std::shared_ptr<ndarray<double>> V1,
      const double* x, double *v) { return false; }

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

inline void particle_tracer::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
  current_t = current_timestep;

  // fprintf(stderr, "#snapshots=%zu\n", this->snapshots.size());
  if (this->snapshots.size() < 2) 
    return; // nothing can be done

  std::shared_ptr<ndarray<double>> V0 = snapshots[0]->get_ptr<double>("vector"),
                                   V1 = snapshots[1]->get_ptr<double>("vector");

  // fprintf(stderr, "#particles=%zu\n", this->particles.size());
  this->parallel_for_container(trajectories, [&](feature_curve_set_t::iterator it) {
  // for (auto &kv : trajectories) {
    // auto &traj = kv.second;
    auto &traj = it->second;
    const auto last_point = traj.back();

    // auto &p = particles[i];
    // double x[3] = {p.x[0], p.x[1], p.t};
   
    double x[nd_+1], v[nd_+1];
    for (auto k = 0; k < nd_; k ++)
      x[k] = last_point.x[k];
    x[nd_] = last_point.t;

    for (auto k = 0; k < nsteps_per_interval; k ++) {
      if (k % nsteps_per_checkpoint == 0) {
        feature_point_t p;
        for (auto k = 0; k < nd_; k ++) {
          p.x[k] = x[k];
          p.v[k] = v[k];
        }
        p.t = x[nd_];
        traj.push_back(p);
      }
      bool succ = rk4<double>(nd_+1, x, [&](const double *x, double *v) { return eval_v(V0, V1, x, v); }, delta(), v);
      if (!succ) 
        break;
    }
        
    feature_point_t p;
    for (auto k = 0; k < nd_; k ++) {
      p.x[k] = x[k];
      p.v[k] = v[k];
    }
    p.t = current_t + current_delta_t; // x[2];
    traj.push_back(p);

    // fprintf(stderr, "%f, %f, %f\n", p.x[0], p.x[1], p.t);
  }, FTK_THREAD_PTHREAD, get_number_of_threads());
}


} // namespace ftk

#endif
