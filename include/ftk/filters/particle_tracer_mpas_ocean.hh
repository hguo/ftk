#ifndef _FTK_PARTICLE_TRACER_MPAS_OCEAN_HH
#define _FTK_PARTICLE_TRACER_MPAS_OCEAN_HH

#include <ftk/ndarray.hh>
#include <ftk/mesh/lattice_partitioner.hh>
#include <ftk/filters/particle_tracer.hh>
#include <ftk/filters/mpas_ocean_tracker.hh>
#include <ftk/utils/gather.hh>
#include <ftk/numeric/wachspress_interpolation.hh>

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
  // 1. find the cell and number of vertices
  // 2. apply wachpress interpolation
  
  int cell_i = m->locate_cell_i(x);
  assert(cell_i < m->n_cells());

  static const int max_nverts = 10;
  int verts_i[max_nverts];
  const int nverts = m->verts_i_on_cell_i(cell_i, verts_i);

  double Xv[max_nverts][3]; // coordinates of vertices
  m->verts_i_coords(nverts, verts_i, Xv);

  double Vv[max_nverts][3]; // velocities on vertices
  for (int i = 0; i < nverts; i ++)
    for (int k = 0; k < 3; k ++)
      Vv[i][k] = V->at(k, i);

  double omega[max_nverts] = {0};
  wachspress_weights(nverts, Xv, x, omega);

  // fprintf(stderr, "omega=%f, %f, %f, %f, %f, %f, %f\n", 
  //     omega[0], omega[1], omega[2], omega[3], omega[4], omega[5], omega[6]);
 
  for (int i = 0; i < nverts; i ++)
    v[i] = omega[i] * V->at(i, verts_i[i]);

  return true;

#if 0
  for (size_t k = 0; k < 3; k ++)
    v[k] = V->get(k, 0, cellId) * 1e5;
  
  fprintf(stderr, "x=%f, %f, %f, %f, cellId=%zu, v=%f, %f, %f\n", 
      x[0], x[1], x[2], x[3], cellId, v[0], v[1], v[2]);
#endif
}

}

#endif
