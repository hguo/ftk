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
  int nverts, nedges, ncells, nlayers, nattrs, max_edges;

  // mesh
  double *d_Xc, *d_Xe, *d_Xv; // cell/edge/vertex coordinates
  int *d_nedges_on_cell, 
      *d_cells_on_cell,
      *d_cells_on_edge,
      *d_cells_on_vert,
      *d_edges_on_cell,
      *d_verts_on_cell;

  // c2v interpolants
  double *d_c2v_interpolants;
  bool *d_vert_on_boundary;
  double *dcw; // a device buffer for c2v interpolation
  double *dew; // a device buffer for e2c interpolation

  // e2c interpolants
  double *d_e2c_interpolants;

  // time-varying data
  double *d_V[2], *d_Vv[2], *d_zTop[2], *d_A[2]; // velocity, verticalVelocity, zTop, and more
  double **dd_V, **dd_Vv, **dd_zTop, **dd_A; // device pointers to pointers
  double T[2];

  // particle data
  int nparticles;
  ftk::feature_point_lite_t *hparts = NULL, *dparts = NULL;
  ftk::feature_point_lite_t *htrajs = NULL, *dtrajs = NULL;

} mop_ctx_t;

void mop_create_ctx(mop_ctx_t **c_, int device=0);
void mop_destroy_ctx(mop_ctx_t **c_);

void mop_load_mesh(mop_ctx_t *c,
    const int ncells,
    const int nedges,
    const int nverts, 
    // const int nlayers, 
    const int max_edges,
    const int nch,
    const double *Xc,
    const double *Xv,
    const int *n_edges_on_cell, 
    const int *cells_on_cell,
    const int *cells_on_edge,
    const int *cells_on_vert,
    const int *edges_on_cell,
    const int *verts_on_cell);

void mop_set_nlayers(mop_ctx_t *c, int n) {c->nlayers = n;}

void mop_load_e2c_interpolants(mop_ctx_t *c,
    const double *p);

void mop_load_data(mop_ctx_t *c, 
    const double *V, 
    const double *Vv, 
    const double *zTop,
    const double *A);

void mop_load_data_with_normal_velocity(mop_ctx_t *c,
    const double t,
    const double *V, // normal velocity
    const double *Vv, 
    const double *zTop,
    const double *A[]);

void mop_load_data_cw(mop_ctx_t *c,
    const double t, // time
    const double *V, 
    const double *Vv, 
    const double *zTop,
    const double *A);

void mop_load_particles(mop_ctx_t *c, 
    const int n,
    ftk::feature_point_lite_t *);

void mop_execute(mop_ctx_t *c, 
    const double T,
    const int nsteps, 
    const int nsubsteps);

#endif
