#ifndef _FTK_PARTICLE_TRACER_2D_REGULAR_HH
#define _FTK_PARTICLE_TRACER_2D_REGULAR_HH

#include <ftk/config.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/clamp.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
// #include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/adjugate.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/geometry/curve2vtk.hh>
#include <ftk/filters/particle_tracer_regular.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/utils/gather.hh>

namespace ftk {

struct particle_tracer_2d_regular : public particle_tracer_regular
{
  particle_tracer_2d_regular(diy::mpi::communicator comm) : 
    particle_tracer_regular(comm, 2), 
    tracker(comm)
  {}
  virtual ~particle_tracer_2d_regular() {}

  void initialize_particles_at_grid_points();
  void update_timestep();

protected:
  bool eval_v(std::shared_ptr<ndarray<double>> V0,
      std::shared_ptr<ndarray<double>> V1,
      const double* x, double *v);
};

inline bool particle_tracer_2d_regular::eval_v(
    std::shared_ptr<ndarray<double>> V0,
    std::shared_ptr<ndarray<double>> V1,
    const double *x, double *v)
{
  if (V0 && V1) { // double time step
    const double t = (x[2] - current_t) / current_delta_t;
    if (t < 0.0 || t > 1.0) return false; // out of temporal bound

    const double w0 = (1.0 - t), w1 = t;
    double v0[2], v1[2];

    const bool b0 = V0->mlerp(x, v0); // current
    const bool b1 = V1->mlerp(x, v1);

    if (b0 && b1) {
      v[0] = w0 * v0[0] + w1 * v1[0];
      v[1] = w0 * v0[1] + w1 * v1[1];
      v[2] = 1.0;
      
      // fprintf(stderr, "x=%f, %f, %f, t=%f, v=%f, %f, %f\n", x[0], x[1], x[2], t, v[0], v[1], v[2]);
      return true;
    } else 
      return false;
  } else if (V0) { // single time step
    return V0->mlerp(x, v);
  } else // no timestep available
    return false;
}

inline void particle_tracer_2d_regular::initialize_particles_at_grid_points()
{
#if 0 // specific points
  feature_curve_t curve;
  feature_point_t p;
  p.x[0] = 1.0; 
  p.x[1] = 1.0;
  // particles.push_back(p);
  curve.push_back(p);

  trajectories.add(curve);
#endif

#if 1
  this->m.element_for(0,
      [&](simplicial_regular_mesh_element e) {
        feature_curve_t curve;
        feature_point_t p;
        p.x[0] = e.corner[0];
        p.x[1] = e.corner[1];
        curve.push_back(p);
        trajectories.add(curve);
      }, 
      FTK_XL_NONE, 1);

  fprintf(stderr, "#trajectories=%zu\n", trajectories.size());
#endif
}

inline void particle_tracer_2d_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
  current_t = current_timestep;

  // fprintf(stderr, "#snapshots=%zu\n", this->snapshots.size());
  if (this->snapshots.size() < 2) 
    return; // nothing can be done

  std::shared_ptr<ndarray<double>> V0 = snapshots[0]->get_ptr<double>("vector"),
                                   V1 = snapshots[1]->get_ptr<double>("vector");

  const int nsteps = 512;
  const int nsteps_per_checkpoint = 128;
  const double delta = 1.0 / nsteps;

  // fprintf(stderr, "#particles=%zu\n", this->particles.size());
  for (auto &kv : trajectories) {
    auto &traj = kv.second;
    const auto last_point = traj.back();

    // auto &p = particles[i];
    // double x[3] = {p.x[0], p.x[1], p.t};
    double x[3] = {last_point.x[0], last_point.x[1], last_point.t};
    double v[3];

    for (auto k = 0; k < nsteps; k ++) {
      if (k % nsteps_per_checkpoint == 0) {
        feature_point_t p;
        p.x[0] = x[0];
        p.x[1] = x[1];
        p.t = x[2];
        p.v[0] = v[0];
        p.v[1] = v[1];
        p.v[2] = v[2];
        traj.push_back(p);
      }
      bool succ = rk4<3, double>(x, [&](const double *x, double *v) { return eval_v(V0, V1, x, v); }, delta, v);
      if (!succ) 
        break;
    }

    feature_point_t p;
    p.x[0] = x[0];
    p.x[1] = x[1];
    p.v[0] = v[0];
    p.v[1] = v[1];
    p.v[2] = v[2];
    p.t = current_t + current_delta_t; // x[2];
    traj.push_back(p);

    // fprintf(stderr, "%f, %f, %f\n", p.x[0], p.x[1], p.t);
  }
}

}

#endif
