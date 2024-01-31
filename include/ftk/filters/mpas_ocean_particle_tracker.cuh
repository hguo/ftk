#ifndef _MPAS_OCEAN_PARTICLE_TRACKER_CUH
#define _MPAS_OCEAN_PARTICLE_TRACKER_CUH

#include <vector>
#include <set>
#include <ftk/features/feature_point_lite.hh>
#include <ftk/basic/kd_lite.hh>

typedef struct {
  // gpu
  int device;

  // precision
  bool prec_compute, prec_mesh, prec_var; // true: double prec

  // mesh info
  int nverts, nedges, ncells, nlayers, nattrs, max_edges;

  // kd-lite
  int *d_kdheap;

  // mesh
  void *d_Xc, *d_Xe, *d_Xv; // cell/edge/vertex coordinates
  int *d_nedges_on_cell, 
      *d_cells_on_cell,
      *d_cells_on_edge,
      *d_cells_on_vert,
      *d_edges_on_cell,
      *d_verts_on_cell;

  // c2v interpolants
  void *d_c2v_interpolants;
  bool *d_vert_on_boundary;
  void *dcw; // a device buffer for c2v interpolation
  void *dew; // a device buffer for e2c interpolation

  // e2c interpolants
  void *d_e2c_interpolants;

  // time-varying data
  void *d_V[2], *d_Vv[2], *d_Z[2], *d_A[2]; // velocity, verticalVelocity, Z, and more
  void **dd_V, **dd_Vv, **dd_Z, **dd_A; // device pointers to pointers
  double T[2]; // time

  // particle data
  int nparticles;
  ftk::feature_point_lite_t *hparts = NULL, *dparts = NULL;
  ftk::feature_point_lite_t *htrajs = NULL, *dtrajs = NULL;

} mop_ctx_t;

void mop_create_ctx(mop_ctx_t **c_, int device, bool prec_compute, bool prec_mesh, bool prec_var);
void mop_destroy_ctx(mop_ctx_t **c_);

void mop_load_mesh(mop_ctx_t *c,
    const int ncells,
    const int nedges,
    const int nverts, 
    // const int nlayers, 
    const int max_edges,
    const int nch,
    const void *Xc,
    const void *Xv,
    const int *n_edges_on_cell, 
    const int *cells_on_cell,
    const int *cells_on_edge,
    const int *cells_on_vert,
    const int *edges_on_cell,
    const int *verts_on_cell);

void mop_set_nlayers(mop_ctx_t *c, int n) {c->nlayers = n;}

void mop_load_kdheap(mop_ctx_t *c, int *heap);

void mop_load_e2c_interpolants(mop_ctx_t *c,
    const void *p);

void mop_seed_latlonz(mop_ctx_t *c,
      const int nlat, const double lat0, const double lat1,
      const int nlon, const double lon0, const double lon1,
      const int nz,   const double z0,   const double z1);

void mop_load_data(mop_ctx_t *c, 
    const void *V, 
    const void *Vv, 
    const void *Z,
    const void *A);

void mop_load_data_with_normal_velocity(mop_ctx_t *c,
    const double t,
    const void *V, // normal velocity
    const void *Vv, 
    const void *z,
    bool isZMid, // true if z is zMid; otherwise z is Z
    const void *layerThickness, // if layerThickness is not null, Z will be derived
    const void *A[]);

void mop_load_data_cw(mop_ctx_t *c,
    const double t, // time
    const void *V, 
    const void *Vv, 
    const void *Z,
    const void *A);

void mop_load_particles(mop_ctx_t *c, 
    const int n,
    ftk::feature_point_lite_t *);

void mop_execute(mop_ctx_t *c, 
    const double T,
    const int nsteps, 
    const int nsubsteps);

#endif
