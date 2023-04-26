#ifndef _FTK_PARTICLE_TRACER_MPAS_OCEAN_HH
#define _FTK_PARTICLE_TRACER_MPAS_OCEAN_HH

#include <ftk/ndarray.hh>
#include <ftk/mesh/lattice_partitioner.hh>
#include <ftk/filters/particle_tracer.hh>
#include <ftk/filters/mpas_ocean_tracker.hh>
#include <ftk/utils/gather.hh>

namespace ftk {

struct particle_tracer_mpas_ocean : public particle_tracer, public mpas_ocean_tracker
{
  particle_tracer_mpas_ocean(diy::mpi::communicator comm, 
      std::shared_ptr<simplicial_mpas_2d_mesh<>> m) : 
    particle_tracer(comm, 3), 
    mpas_ocean_tracker(comm, m), 
    tracker(comm) 
  {} 

  virtual ~particle_tracer_mpas_ocean() {}

  void initialize_particles_at_grid_points();

protected:
  bool eval_v(std::shared_ptr<ndarray<double>> V0,
      std::shared_ptr<ndarray<double>> V1,
      const double* x, double *v);
};

////
inline void particle_tracer_mpas_ocean::initialize_particles_at_grid_points()
{
#if 0
  this->m.element_for_ordinal(0, 0,
      [&](simplicial_mpas_ocean_mesh_element e) {
        feature_curve_t curve;
        feature_point_t p;
        for (auto k = 0; k < std::min(size_t(3), e.corner.size()-1); k ++) { // m is (n+1)-d mesh
          p.x[k] = e.corner[k];  
        }
        // fprintf(stderr, "%f, %f, %f\n", p.x[0], p.x[1], p.x[2]);
        curve.push_back(p);
        trajectories.add(curve);
      }, 
      FTK_XL_NONE, 1);
#endif
  fprintf(stderr, "#trajectories=%zu\n", trajectories.size());
}

inline bool particle_tracer_mpas_ocean::eval_v(
    std::shared_ptr<ndarray<double>> V0,
    std::shared_ptr<ndarray<double>> V1,
    const double *x, double *v)
{
#if 0
  if (V0 && V1) { // double time step
    const double t = (x[nd_] - current_t) / current_delta_t;
    if (t < 0.0 || t > 1.0) return false; // out of temporal bound

    const double w0 = (1.0 - t), w1 = t;
    double v0[4], v1[4]; // up to 4d for now

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
#endif
  return false;
}

}

#endif
