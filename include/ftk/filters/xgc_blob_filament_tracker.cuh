#ifndef _XGC_BLOB_FILAMENT_TRACKER_CUH
#define _XGC_BLOB_FILAMENT_TRACKER_CUH

#include <vector>
#include <set>
#include <ftk/numeric/xgc_interpolant.hh>
#include <ftk/features/feature_point_lite.hh>
#include <ftk/mesh/bvh2d.hh>

typedef struct {
  int m2n0, m2n1, m2n2, max_vertex_triangles;
  int nphi = 16, iphi = 1, vphi = 1;

  // mesh
  double *d_m2coords, *d_m2invdet = NULL;
  int *d_m2edges, *d_m2tris;
  int *d_vertex_triangles = NULL;
  double *d_psin; // normalized psi

  // bvh
  bvh2d_node_t<> *d_bvh = NULL;

  // interpolants
  ftk::xgc_interpolant_t<> *d_interpolants = NULL;

  // smoothing kernel
  double sigma;
  int *d_kernel_nodes;
  double *d_kernel_values;
  size_t *d_kernel_lengths, *d_kernel_offsets;

  double *d_scalar_in;
  double *d_scalar[2], *d_vector[2], *d_jacobian[2];

  ftk::feature_point_lite_t *hcps = NULL, *dcps = NULL;
  unsigned long long hncps = 0, *dncps = NULL;
  size_t bufsize; //  = 512 * 1024 * 1024; // 512 MB of buffer
  int device;

  double factor; // scaling factor
  
  // for poincare plot
  bool poincare = false;
  double *d_apars = NULL, *d_apars_upsample = NULL, *h_apars_upsample; // apars (nphi) and its upsampled (nphi*iphi) version
  double *d_gradAs = NULL, *d_gradAs_cw = NULL; // 2D gradient of upsampled apars, vertexwise and cellwise
  double *d_bfield = NULL, *d_bfield0 = NULL, *d_curl_bfield0 = NULL;
  double **ddp_deltaB = NULL, **hdp_deltaB;
  // double *d_deltaB = NULL, *h_deltaB; // (upsampled) deltaB
  // double *d_seeds = NULL;
  double *d_poincare_psin, *h_poincare_psin;
  int nseeds, nsteps;

  bool retrieve_apars_upsample;
} xft_ctx_t;

void xft_create_ctx(xft_ctx_t **c_, int device=0, int buffer_size_in_mb=512);
void xft_create_poincare_ctx(xft_ctx_t **c_, int nseeds, int nsteps, int device=0);
void xft_destroy_ctx(xft_ctx_t **c_);
void xft_load_mesh(xft_ctx_t *c,
    int nphi, int iphi, int vphi,
    int m2n0, int m2n1, int m2n2,
    const double *m2coords, const int *m2edges, const int *m2tris);
void xft_load_vertex_triangles(xft_ctx_t *c, const std::vector<std::set<int>>& vertex_triangles);
void xft_load_bvh(xft_ctx_t *c, const std::vector<bvh2d_node_t<int, double>>& bvh);
void xft_derive_interpolants(xft_ctx_t *c);
void xft_load_interpolants(xft_ctx_t *c, const std::vector<std::vector<ftk::xgc_interpolant_t<>>> &interpolants);
void xft_load_smoothing_kernel(xft_ctx_t *c, double sigma, const std::vector<std::vector<std::tuple<int, double>>>& kernels);
void xft_load_psin(xft_ctx_t *c, const double *psin);
void xft_execute(xft_ctx_t *c, int scope, int current_timestep);
void xft_load_scalar_data(xft_ctx_t *c, const double *scalar);
void xft_load_data(xft_ctx_t *c, 
    const double *scalar, 
    const double *vector, 
    const double *jacobian);
void xft_swap(xft_ctx_t *c);

// poincare
void xft_load_magnetic_field(xft_ctx_t *c, 
    const double *bfield, 
    const double *bfield0, 
    const double *curl_bfield0);

void xft_load_apars(xft_ctx_t *c, 
    const double *apars);

void xft_load_apars_upsample(xft_ctx_t *c, 
    const double *apars_upsample);

void xft_retrieve_deltaB(xft_ctx_t *c);

void xft_compute_poincare_plot(xft_ctx_t *c,
    const double *seeds, bool use_static_b = false, 
    int dir = 1);

void xft_compute_poincare_psin(xft_ctx_t *c);

#endif
