#ifndef _MPAS_OCEAN_PARTICLE_TRACKER_CUH
#define _MPAS_OCEAN_PARTICLE_TRACKER_CUH

#include <vector>
#include <set>
#include <ftk/features/feature_point_lite.hh>
#include <ftk/mesh/bvh2d.hh>

typedef struct {
  // gpu
  int device;

  // mesh info
  int ncells, nlayers, nverts, max_edges, nattrs;

  // mesh
  double *d_Xc, *d_Xv; // cell/vertex coordinates
  int *d_nedges_on_cell, 
      *d_cells_on_cell,
      *d_verts_on_cell;

  // time-varying data
  double *d_V[2], // velocity fields of adjacent timesteps
         *d_Vv[2], // vertical velocities
         *d_zTop[2], // top layer depth
         *d_A[2]; // scalar attrs

  // particle data
  int nparticles;
  ftk::feature_point_lite_t *hcps = NULL, *dcps = NULL;

} mop_ctx_t;

void mop_create_ctx(mop_ctx_t **c_, int device=0);
void mop_destroy_ctx(mop_ctx_t **c_);
void mop_load_mesh(mop_ctx_t *c,
    const int ncells, 
    const int nlayers, 
    const int nverts, 
    const int max_edges,
    const int nattrs,
    const double *n_edges_on_cell, 
    const double *cells_on_cell,
    const double *verts_on_cell);

void mop_load_data(mop_ctx_t *c,
    const double *V, 
    const double *Vv, 
    const double *zTop, 
    const double *A);

void mop_execute(mop_ctx_t *c, int scope, int current_timestep);
void mop_swap(mop_ctx_t *c);

#endif
