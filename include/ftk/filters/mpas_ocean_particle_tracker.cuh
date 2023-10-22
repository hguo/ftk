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
  int ncells, nlayers, nverts, nattrs, max_edges;

  // mesh
  double *d_Xc, *d_Xv; // cell/vertex coordinates
  int *d_nedges_on_cell, 
      *d_cells_on_cell,
      *d_verts_on_cell,
      *d_cells_on_vert;

  // c2v interpolants
  double *d_c2v_interpolants;
  bool *d_cell_on_boundary;

  // time-varying data
  double *d_V[2], *d_Vv[2], *d_zTop[2], *d_A[2]; // velocity, verticalVelocity, zTop, and more

  // particle data
  int nparticles;
  ftk::feature_point_lite_t *hparts = NULL, *dparts = NULL;
  ftk::feature_point_lite_t *htrajs = NULL, *dtrajs = NULL;

} mop_ctx_t;

void mop_create_ctx(mop_ctx_t **c_, int device=0);
void mop_destroy_ctx(mop_ctx_t **c_);

void mop_load_mesh(mop_ctx_t *c,
    const int ncells, 
    const int nlayers, 
    const int nverts, 
    const int max_edges,
    const int nch,
    const double *Xc,
    const double *Xv,
    const int *n_edges_on_cell, 
    const int *cells_on_cell,
    const int *verts_on_cell,
    const int *cells_on_vert);

void mop_load_data(mop_ctx_t *c, 
    const double *V, 
    const double *Vv, 
    const double *zTop,
    const double *A);

void mop_load_particles(mop_ctx_t *c, 
    const int n,
    ftk::feature_point_lite_t *);

void mop_execute(mop_ctx_t *c, int current_timestep);
void mop_swap(mop_ctx_t *c);

#endif
