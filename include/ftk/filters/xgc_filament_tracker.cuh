#ifndef _XGC_FILAMENT_TRACKER_CUH
#define _XGC_FILAMENT_TRACKER_CUH

#include <vector>
#include <ftk/numeric/xgc_interpolant.hh>
#include <ftk/filters/feature_point_lite.hh>

typedef struct {
  int m2n0, m2n1, m2n2;
  int nphi = 16, iphi = 1, vphi = 1;

  double *d_m2coords;
  int *d_m2edges, *d_m2tris;
  ftk::xgc_interpolant_t<> **d_ptr_interpolants = NULL, *d_interpolants = NULL;

  double *d_scalar[2] = {NULL, NULL}, *d_vector[2] = {NULL, NULL}, *d_jacobian[2] = {NULL, NULL};

  ftk::feature_point_lite_t *hcps = NULL, *dcps = NULL;
  unsigned long long hncps = 0, *dncps = NULL;
  const size_t bufsize = 512 * 1024 * 1024; // 512 MB of buffer
} xft_ctx_t;

void xft_create_ctx(xft_ctx_t **c_);
void xft_destroy_ctx(xft_ctx_t **c_);
void xft_load_mesh(xft_ctx_t *c,
    int nphi, int iphi, int vphi,
    int m2n0, int m2n1, int m2n2,
    const double *m2coords, const int *m2edges, const int *m2tris);
void xft_load_interpolants(xft_ctx_t *c, const std::vector<std::vector<ftk::xgc_interpolant_t<>>> &interpolants);
void xft_execute(xft_ctx_t *c, int scope, int current_timestep);
void xft_load_data(xft_ctx_t *c, 
    const double *scalar, 
    const double *vector, 
    const double *jacobian);


#endif
