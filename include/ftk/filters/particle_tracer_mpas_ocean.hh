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
      std::shared_ptr<mpas_mesh<>> m) : 
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
  for (auto i = 0; i < m->n_cells(); i += 1000000) {
        feature_curve_t curve;
        feature_point_t p;
        for (auto k = 0; k < 3; k ++)
          p.x[k] = m->xyzCells(i, k); 
        
        curve.push_back(p);
        trajectories.add(curve);
  }
  
  fprintf(stderr, "#trajectories=%zu\n", trajectories.size());
}

inline bool particle_tracer_mpas_ocean::eval_v(
    std::shared_ptr<ndarray<double>> V,
    const double *x, double *v)
{
  size_t cellId = m->locate_cell(x);
  assert(cellId < m->n_cells());

  for (size_t k = 0; k < 3; k ++)
    v[k] = V->get(k, 0, cellId) * 1e5;
  
  fprintf(stderr, "x=%f, %f, %f, %f, cellId=%zu, v=%f, %f, %f\n", 
      x[0], x[1], x[2], x[3], cellId, v[0], v[1], v[2]);

  return true;
  // return false;
  
  // 1. find the nearest vertex
  // 2. barycentric interpolation
}

}

#endif
