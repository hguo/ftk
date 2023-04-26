#ifndef _FTK_PARTICLE_TRACER_REGULAR_HH
#define _FTK_PARTICLE_TRACER_REGULAR_HH

#include <ftk/ndarray.hh>
#include <ftk/mesh/lattice_partitioner.hh>
#include <ftk/filters/particle_tracer.hh>
#include <ftk/filters/regular_tracker.hh>
#include <ftk/utils/gather.hh>

namespace ftk {

struct particle_tracer_regular : public particle_tracer, public regular_tracker
{
  particle_tracer_regular(diy::mpi::communicator comm, int nd) : 
    particle_tracer(comm, nd), 
    regular_tracker(comm, nd), 
    tracker(comm) 
  {} 

  virtual ~particle_tracer_regular() {}

  void initialize_particles_at_grid_points();
  void update_timestep();

protected:
  bool eval_v(std::shared_ptr<ndarray<double>> V0,
      std::shared_ptr<ndarray<double>> V1,
      const double* x, double *v);
};

////
inline void particle_tracer_regular::initialize_particles_at_grid_points()
{
  this->m.element_for_ordinal(0, 0,
      [&](simplicial_regular_mesh_element e) {
        feature_curve_t curve;
        feature_point_t p;
        for (auto k = 0; k < std::min(size_t(3), e.corner.size()-1); k ++) { // m is (n+1)-d mesh
          p.x[k] = e.corner[k];  
        }
        curve.push_back(p);
        trajectories.add(curve);
      }, 
      FTK_XL_NONE, 1);

  fprintf(stderr, "#trajectories=%zu\n", trajectories.size());
}

inline bool particle_tracer_regular::eval_v(
    std::shared_ptr<ndarray<double>> V0,
    std::shared_ptr<ndarray<double>> V1,
    const double *x, double *v)
{
  if (V0 && V1) { // double time step
    const double t = (x[nd_] - current_t) / current_delta_t;
    if (t < 0.0 || t > 1.0) return false; // out of temporal bound

    const double w0 = (1.0 - t), w1 = t;
    double v0[2], v1[2];

    const bool b0 = V0->mlerp(x, v0); // current
    const bool b1 = V1->mlerp(x, v1);

    if (b0 && b1) {
      for (auto k = 0; k < nd(); k ++)
        v[k] = w0 * v0[k] + w1 * v1[k];
      v[nd()] = 1.0; // time
      
      // fprintf(stderr, "x=%f, %f, %f, t=%f, v=%f, %f, %f\n", x[0], x[1], x[2], t, v[0], v[1], v[2]);
      return true;
    } else 
      return false;
  } else if (V0) { // single time step
    return V0->mlerp(x, v);
  } else // no timestep available
    return false;
}

inline void particle_tracer_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
  current_t = current_timestep;

  // fprintf(stderr, "#snapshots=%zu\n", this->snapshots.size());
  if (this->snapshots.size() < 2) 
    return; // nothing can be done

  std::shared_ptr<ndarray<double>> V0 = snapshots[0]->get_ptr<double>("vector"),
                                   V1 = snapshots[1]->get_ptr<double>("vector");

  // fprintf(stderr, "#particles=%zu\n", this->particles.size());
  // for (auto &kv : trajectories) {
  //   auto &traj = kv.second;
  this->parallel_for_container(trajectories, [&](feature_curve_set_t::iterator it) {
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
  }, FTK_THREAD_PTHREAD, 1);
}


}

#endif
