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
  bool eval_v(std::shared_ptr<ndarray<double>> V,
      const double* x, double *v);
};

////
inline void particle_tracer_mpas_ocean::initialize_particles_at_grid_points()
{
  m->element_for(0, [&](int i) {
        double x[3];
        m->get_coords(i, x);

        feature_curve_t curve;
        feature_point_t p;
        for (auto k = 0; k < 3; k ++)
          p.x[k] = x[k]; 
        
        curve.push_back(p);
        trajectories.add(curve);
  });  // FTK_XL_NONE, 1);
  fprintf(stderr, "#trajectories=%zu\n", trajectories.size());
}

inline bool particle_tracer_mpas_ocean::eval_v(
    std::shared_ptr<ndarray<double>> V,
    const double *x, double *v)
{
  return false;
  
  // 1. find the nearest vertex
  // 2. barycentric interpolation
}

}

#endif
