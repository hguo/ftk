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
  bool eval_v(const double* x, double *v);
};

inline bool particle_tracer_2d_regular::eval_v(
    const double *x, double *v)
{
  const double t = (x[3] - current_t) / current_delta_t;
  const double w0 = (1.0 - t), w1 = t;
  double v0[2], v1[2];

  field_data_snapshots[0].vector.mlerp(x, v0); // current
  field_data_snapshots[1].vector.mlerp(x, v1);

  v[0] = w0 * v0[0] + w1 * v1[0];
  v[1] = w0 * v0[1] + w1 * v1[1];

  return false;
}

inline void particle_tracer_2d_regular::initialize_particles_at_grid_points()
{
  this->m.element_for(0,
      [&](simplicial_regular_mesh_element e) {
        feature_point_lite_t p;
        p.x[0] = e.corner[0];
        p.x[1] = e.corner[1];
        particles.push_back(p);
      }, 
      FTK_XL_NONE, 1);
}

inline void particle_tracer_2d_regular::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
}

}

#endif
