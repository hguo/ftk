#ifndef _XGC_BLOB_FILAMENT_TRACKER_CUH
#define _XGC_BLOB_FILAMENT_TRACKER_CUH

#include <vector>
#include <set>
#include <ftk/numeric/xgc_interpolant.hh>
#include <ftk/features/feature_point_lite.hh>
#include <ftk/mesh/bvh2d.hh>

static const int maxgpus = 8;

typedef struct {
  int m2n0, m2n1, m2n2, max_vertex_triangles;
  int nphi = 16, iphi = 1, vphi = 1;

  // mesh
  double *d_m2coords[maxgpus], *d_m2invdet[maxgpus];
  int *d_m2edges[maxgpus], *d_m2tris[maxgpus];
  int *d_vertex_triangles[maxgpus];
  double *d_psin[maxgpus]; // normalized psi

  // bvh
  bvh2d_node_t<> *d_bvh[maxgpus];

  // interpolants
  ftk::xgc_interpolant_t<> *d_interpolants[maxgpus];

  // smoothing kernel
  double sigma;
  int *d_kernel_nodes[maxgpus];
  double *d_kernel_values[maxgpus];
  size_t *d_kernel_lengths[maxgpus], *d_kernel_offsets[maxgpus];

  double *d_scalar_in[maxgpus];
  double *d_scalar[maxgpus][2], *d_vector[maxgpus][2], *d_jacobian[maxgpus][2];

  ftk::feature_point_lite_t *hcps, *dcps[maxgpus];
  unsigned long long hncps = 0, *dncps[maxgpus];
  size_t bufsize; //  = 512 * 1024 * 1024; // 512 MB of buffer
  int ndevices;

  double factor; // scaling factor
  
  // for poincare plot
  bool poincare = false;
  double *d_apars[maxgpus], *d_apars_upsample[maxgpus]; // apars (nphi) and its upsampled (nphi*iphi) version
  double *d_gradAs[maxgpus], *d_gradAs_cw[maxgpus]; // 2D gradient of upsampled apars, vertexwise and cellwise
  double *d_bfield[maxgpus], *d_bfield0[maxgpus], *d_curl_bfield0[maxgpus];
  double *d_deltaB[maxgpus]; // (upsampled) deltaB
  double *d_seeds[maxgpus];
  double *d_poincare_psin[maxgpus], *h_poincare_psin;
  int nseeds, nsteps;
} xft_ctx_t;

void xft_create_ctx(xft_ctx_t **c_, int buffer_size_in_mb=512);
void xft_create_poincare_ctx(xft_ctx_t **c_, int nseeds, int nsteps);
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

void xft_compute_poincare_plot(xft_ctx_t *c,
    const double *seeds, bool use_static_b = false, 
    int dir = 1);

void xft_compute_poincare_psin(xft_ctx_t *c);

#endif
