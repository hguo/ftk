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

protected:
  virtual bool eval_v(std::shared_ptr<ndarray<double>> V, // single timestep vector field
      const double *x, double *v);
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
        // fprintf(stderr, "%f, %f, %f\n", p.x[0], p.x[1], p.x[2]);
        curve.push_back(p);
        trajectories.add(curve);
      }, 
      FTK_XL_NONE, 1);

  fprintf(stderr, "#trajectories=%zu\n", trajectories.size());
}
  
bool particle_tracer_regular::eval_v(std::shared_ptr<ndarray<double>> V, // single timestep vector field
    const double *x, double *v)
{
  return V->mlerp(x, v);
}

}

#endif
