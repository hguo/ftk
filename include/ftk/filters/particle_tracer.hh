#ifndef _FTK_PARTICLE_TRACER_HH
#define _FTK_PARTICLE_TRACER_HH

#include <ftk/config.hh>
#include <ftk/filters/tracker.hh>

namespace ftk {

enum {
  PARTICLE_TRACER_INTEGRATOR_RK1 = 0,
  PARTICLE_TRACER_INTEGRATOR_RK4 = 1,
  PARTICLE_TRACER_INTEGRATOR_SPHERICAL_RK1 = 5,
  PARTICLE_TRACER_INTEGRATOR_SPHERICAL_RK1_WITH_VERTICAL_VELOCITY = 6,
  PARTICLE_TRACER_INTEGRATOR_SPHERICAL_RK4_WITH_VERTICAL_VELOCITY = 7,
};

struct particle_tracer : public virtual tracker
{
  particle_tracer(diy::mpi::communicator comm, int nd /*spatial dim*/) : tracker(comm), nd_(nd) {}
  virtual ~particle_tracer() {}

  virtual void initialize_particles(size_t n, const double *p, size_t stride=3, bool has_time=false);
  virtual void initialize_particles_at_grid_points(std::vector<int> strides = {}) {};

  void update() {}
  void finalize() {}
  void update_timestep();
  bool advance_timestep();

  virtual void prepare_timestep();

  void write_trajectories(const std::string& filename);
  void write_geo_trajectories(const std::string& filename);

  void set_delta_t(double d) { current_delta_t = d; } // overriding delta t
  void set_nsteps_per_interval(int n) { nsteps_per_interval = n; }
  void set_nsteps_per_checkpoint(int n) { nsteps_per_checkpoint = 16; }

protected:
  virtual bool eval_vt(const double *x, double *v, int *hint = NULL);  // w/ temporal interpolation
  virtual bool eval_v(int t, const double *x, double *v, int *hint = NULL) { return false; } // single timestep

  int nd() const { return nd_; }
  double delta() const { return current_delta_t / nsteps_per_interval; }

  virtual int nch() const {return nd_ + 1;} // number of channels, e.g., 2D vel w/ time will be 3
  virtual std::vector<std::string> scalar_names() const { return {}; }

protected:
  std::shared_ptr<ndarray<double>> V[2];

protected:
  std::vector<feature_point_lite_t> particles;
  feature_curve_set_t trajectories;

  double current_t = 0.0, current_delta_t = 1.0;

  int nsteps_per_interval = 128;
  int nsteps_per_checkpoint = 16;

  const int nd_;
  int integrator = PARTICLE_TRACER_INTEGRATOR_RK4;
};

//// 
inline bool particle_tracer::advance_timestep()
{
  update_timestep();
  pop_field_data_snapshot();

  current_timestep ++;
  return snapshots.size() > 0;
}

inline void particle_tracer::write_geo_trajectories(const std::string& filename)
{
#if FTK_HAVE_VTK
  if (comm.rank() == 0) {
    auto poly = trajectories.to_geo().to_vtp( scalar_names() );
    write_polydata(filename, poly);
  }
#endif
}

inline void particle_tracer::write_trajectories(const std::string& filename)
{
#if FTK_HAVE_VTK
  if (comm.rank() == 0) {
    auto poly = trajectories.to_vtp( scalar_names() );
    write_polydata(filename, poly);
  }
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

inline void particle_tracer::prepare_timestep()
{
  V[0] = snapshots[0]->get_ptr<double>("vector");
  V[1] = snapshots.size() > 1 ? snapshots[1]->get_ptr<double>("vector") : nullptr;
}

inline void particle_tracer::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
  current_t = current_timestep;
  prepare_timestep();

  bool streamlines = false;
  if (this->ntimesteps == 1)
    streamlines = true;

  // fprintf(stderr, "#snapshots=%zu\n", this->snapshots.size());
  if (!streamlines && this->snapshots.size() < 2) 
    return; // nothing can be done

  // fprintf(stderr, "#particles=%zu\n", this->particles.size());
  this->parallel_for_container(trajectories, [&](feature_curve_set_t::iterator it) {
  // for (auto &kv : trajectories) {
    // auto &traj = kv.second;
    auto &traj = it->second;
    const auto last_point = traj.back();

    // auto &p = particles[i];
    // double x[3] = {p.x[0], p.x[1], p.t};
   
    double x[10], v[10]; // some large number for more than just 4 components
    for (auto k = 0; k < nd_; k ++)
      x[k] = last_point.x[k];
    x[nd_] = last_point.t;

    thread_local int hint[2] = {-1, -1};

    for (auto k = 0; k < nsteps_per_interval; k ++) {
      if (k % nsteps_per_checkpoint == 0) {
        feature_point_t p;
        for (auto k = 0; k < nd_; k ++) {
          p.x[k] = x[k];
          p.v[k] = v[k];
        }
        p.t = x[nd_]; //  + delta() * k;
        for (auto k = 0; k < nch() - nd() - 1; k ++)
          p.scalar[k] = v[k + nd() + 1];
        // fprintf(stderr, "v=%f, %f, %f, %f, %f, %f, %f\n", v[0], v[1], v[2], v[3], v[4], v[5], v[6]);
        // p.print(std::cerr, scalar_names()) << std::endl;
        p.id = traj.id;
        traj.push_back(p);
      }
      bool succ = false;

      // fprintf(stderr, "x=%f, %f, t=%f\n", x[0], x[1], x[2]);
      if (integrator == PARTICLE_TRACER_INTEGRATOR_RK4)
        succ = rk4<double>(nd_+1, x, [&](const double *x, double *v) { return eval_vt(x, v, hint); }, delta(), v);
      else if (integrator == PARTICLE_TRACER_INTEGRATOR_SPHERICAL_RK1)
        succ = spherical_rk1<double>(x, [&](const double *x, double *v) { return eval_vt(x, v, hint); }, delta(), v);
      else if (integrator == PARTICLE_TRACER_INTEGRATOR_SPHERICAL_RK1_WITH_VERTICAL_VELOCITY) 
        succ = spherical_rk1_with_vertical_velocity<double>(x, [&](const double *x, double *v) { return eval_vt(x, v, hint); }, delta(), v);
      else if (integrator == PARTICLE_TRACER_INTEGRATOR_SPHERICAL_RK4_WITH_VERTICAL_VELOCITY) 
        succ = spherical_rk4_with_vertical_velocity<double>(x, [&](const double *x, double *v) { return eval_vt(x, v, hint); }, delta(), v);

      if (!succ) 
        break;
    }
        
    feature_point_t p;
    for (auto k = 0; k < nd_; k ++) {
      p.x[k] = x[k];
      p.v[k] = v[k];
    }
    // p.t = current_t + current_delta_t; // x[2];
    p.t = x[nd_];
    for (auto k = 0; k < nch() - nd() - 1; k ++)
      p.scalar[k] = v[k + nd() + 1];
    p.id = traj.id;
    traj.push_back(p);

    // fprintf(stderr, "%f, %f, %f\n", p.x[0], p.x[1], p.t);
  }, FTK_THREAD_PTHREAD, get_number_of_threads());
}

inline bool particle_tracer::eval_vt(
    const double *x, double *v, int *hint)
{
  if (V[0] && V[1]) { // double time step
    const double t = (x[nd_] - current_t) / current_delta_t;
    if (t < 0.0 || t > 1.0) return false; // out of temporal bound

    const double w0 = (1.0 - t), w1 = t;
    double v0[10], v1[10]; // some large number

    const bool b0 = eval_v(0, x, v0, hint);
    const bool b1 = eval_v(1, x, v1, hint);

    if (b0 && b1) {
      for (auto k = 0; k < nch(); k ++)
        v[k] = w0 * v0[k] + w1 * v1[k];
      v[nd()] = 1.0; // time
      
      // fprintf(stderr, "x=%f, %f, t=%f, v=%f, %f, %f\n", x[0], x[1], t, v[0], v[1], v[2]);
      return true;
    } else 
      return false;
  } else if (V[0]) { // single time step
    bool succ = eval_v(0, x, v, hint);
    v[nd()] = 1.0; // time
    return succ;
  } else // no timestep available
    return false;
}

} // namespace ftk

#endif
