#include <nvfunctional>
#include "common.cuh"
#include "mx4.cuh"

// what are needed in device memory:
// - field-following interpolants for virtual poloidal planes
// - triangles, edges, and vertex coordinates in 2D mesh
// - current and next timestep of scalar, vector, and jacobian fields

using namespace ftk;

typedef struct {
  int m2n0, m2n1, m2n2;
  int nphi = 16, iphi = 1, vphi = 1;

  double *d_m2coords;
  int *d_m2edges, *d_m2tris;
  xgc_interpolant_t **d_ptr_interpolants = NULL, *d_interpolants = NULL;

  double *d_scalar[2] = {0}, *d_vector[2] = {0}, *d_jacobian[2] = {0};

  cp_t *hcps = NULL, *dcps = NULL;
  unsigned long long hncps = 0, *dncps = NULL;
  const size_t bufsize = 512 * 1024 * 1024; // 512 MB of buffer
} xft_ctx_t;

typedef xft_ctx_t ctx_t;

template <typename I, typename F>
__device__
bool check_simplex(
    I current_timestep, 
    I i,
    const I nphi, const I iphi, const I vphi,
    const I m2n0, const I m2n1, const I m2n2,
    const F m2coords[], // m2 vertex coordinates
    const I m2edges[], // list of m2 edges
    const I m2tris[], // list of m2 triangles
    const xgc_interpolant_t *const *const interpolants, // interpolants
    const F *const scalar[2], // current and next scalar
    const F *const vector[2], // current and next grad
    const F *const jacobian[2], // current and next jacobian
    cp_t & cp) // WIP: critical points
{
  // typedef ftk::fixed_point<> fp_t;

  const I np = nphi * iphi * vphi;
  const I m3n0 = m2n0 * np;

  I verts[3], t[3], p[3];
  F rzpt[3][4], f[3], v[3][2], j[3][2][2];

  mx4_get_tri(i, verts, np, m2n0, m2n1, m2n2, m2edges, m2tris);

  for (int k = 0; k < 3; k ++) {
    t[k] = verts[k] / m3n0; // time in m4
    const I v3 = verts[k] % m3n0; // vert in m3
    const I v2 = v3 % m2n0; // vert in m2
    p[k] = v3 / m2n0; // poloidal plane

    mx3_get_coords(v3, rzpt[k]);
    rzpt[k][3] = t[k];

    const int iv = (t[k] == current_timestep) ? 0 : 1;
    mx3_interpolate(interpolants, 
        scalar[iv], vector[iv], jacobian[iv], 
        v3, f[k], v[k], j[k]);
  }

  // check if peridocial
  bool b0 = false, b1 = false;
  for (int k = 0; k < 3; k ++) {
    if (p[k] == 0) b0 = true;
    else if (p[k] == np - 1) b1 = true;
  }
  if (b0 && b1) // peridical correction
    for (int k = 0; k < 3; k ++) 
      if (p[k] == 0)
        rzpt[k][2] += np;

  // fp_t vf[3][2];
}

template <typename I, typename F>
__global__
void sweep_simplices(
    int scope, 
    I current_timestep, 
    const I nphi, const I iphi, const I vphi,
    const I m2n0, const I m2n1, const I m2n2,
    const F m2coords[], // m2 vertex coordinates
    const I m2edges[], // list of m2 edges
    const I m2tris[], // list of m2 triangles
    const xgc_interpolant_t *const *const interpolants, // interpolants
    const F * const scalar[2], // current and next scalar
    const F * const vector[2], // current and next grad
    const F * const jacobian[2], // current and next jacobian
    unsigned long long &ncps, 
    cp_t *cps)
{
#if 0
  const I mx3n2 = (3 * m2n0 + 2 * m2n1) * np;

  int tid = getGlobalIdx_3D_1D();
  I i = tid;
  if (scope == scope_interval) i += mx3n2;
  
  cp_t cp;
  bool succ = check_simplex(
      current_timestep, 
      i, 
      nphi, iphi, vphi, 
      m2n0, m2n1, m2n2, 
      m2coords, m2edges, m2tris,
      interpolants, 
      scalar, vector, jacobian, 
      ncps, cps);
  
  if (succ) {
    unsigned long long idx = atomicAdd(&ncps, 1ul);
    cp.tag = tid;
    cps[idx] = cp;
  }
#endif
}

void xft_create_ctx(ctx_t **c_)
{
  *c_ = (ctx_t*)malloc(sizeof(ctx_t));
  ctx_t *c = *c_;
  
  cudaMalloc((void**)&c->dncps, sizeof(unsigned long long));
  cudaMemset(c->dncps, 0, sizeof(unsigned long long));
  checkLastCudaError("[FTK-CUDA] cuda malloc");

  c->hcps = (cp_t*)malloc(c->bufsize);
  cudaMalloc((void**)&c->dcps, c->bufsize);
  checkLastCudaError("[FTK-CUDA] cuda malloc");
}

void xft_destroy_ctx(ctx_t **c_)
{
  ctx_t *c = *c_;

  if (c->d_m2coords != NULL) cudaFree(c->d_m2coords);
  if (c->d_m2edges != NULL) cudaFree(c->d_m2edges);
  if (c->d_m2tris != NULL) cudaFree(c->d_m2tris);

  if (c->d_ptr_interpolants != NULL) cudaFree(c->d_ptr_interpolants);
  if (c->d_interpolants != NULL) cudaFree(c->d_interpolants);

  if (c->d_scalar[0] != NULL) cudaFree(c->d_scalar[0]);
  if (c->d_scalar[1] != NULL) cudaFree(c->d_scalar[1]);
  if (c->d_vector[0] != NULL) cudaFree(c->d_vector[0]);
  if (c->d_vector[1] != NULL) cudaFree(c->d_vector[1]);
  if (c->d_jacobian[0] != NULL) cudaFree(c->d_jacobian[0]);
  if (c->d_jacobian[1] != NULL) cudaFree(c->d_jacobian[1]);
  
  checkLastCudaError("[FTK-CUDA] cuda free");

  free(*c_);
  *c_ = NULL;
}

void xft_execute(ctx_t *c, int scope, int current_timestep)
{
  const int np = c->nphi * c->iphi * c->vphi;
  const int mx3n1 = (2 * c->m2n1 + 2 * c->m2n0) * np;
  const int mx3n2 = (3 * c->m2n0 + 2 * c->m2n1) * np;
  // const int mx4n2 = 3 * mx3n2 + 2 * mx3n1;
  const int mx4n2_ordinal  = mx3n2, 
            mx4n2_interval = 2 * mx3n2 + 2 * mx3n1;

  size_t ntasks;
  if (scope == scope_ordinal) ntasks = mx4n2_ordinal; 
  else ntasks = mx4n2_interval;
  
  const int maxGridDim = 1024;
  const int blockSize = 256;
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;

  if (nBlocks >= maxGridDim) 
    gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else 
    gridSize = dim3(nBlocks);

  sweep_simplices<int, double><<<gridSize, blockSize>>>(
      scope, current_timestep, 
      c->nphi, c->iphi, c->vphi, 
      c->m2n0, c->m2n1, c->m2n2, 
      c->d_m2coords, c->d_m2edges, c->d_m2tris, 
      c->d_ptr_interpolants, 
      c->d_scalar, c->d_vector, c->d_jacobian, 
      *c->dncps, c->dcps);
  checkLastCudaError("[FTK-CUDA] sweep_simplicies");

  cudaMemcpy(&c->hncps, c->dncps, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
  cudaMemcpy(c->hcps, c->dcps, sizeof(cp_t) * c->hncps, cudaMemcpyDeviceToHost);
  
  checkLastCudaError("[FTK-CUDA] cuda memcpy device to host");
}

void xft_load_data(ctx_t *c, 
    const double *scalar, const double *vector, const double *jacobian)
{
  double *dd_scalar;
  if (c->d_scalar[0] == NULL) {
    cudaMalloc((void**)&c->d_scalar[0], sizeof(double) * c->m2n0);
    dd_scalar = c->d_scalar[0];
  } else if (c->d_scalar[1] == NULL) {
    cudaMalloc((void**)&c->d_scalar[1], sizeof(double) * c->m2n0);
    dd_scalar = c->d_scalar[1];
  } else {
    std::swap(c->d_scalar[0], c->d_scalar[1]);
    dd_scalar = c->d_scalar[1];
  }
  cudaMemcpy(dd_scalar, scalar, sizeof(double) * c->m2n0, 
      cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] loading scalar field data");
 
  /// 
  double *dd_vector;
  if (c->d_vector[0] == NULL) {
    cudaMalloc((void**)&c->d_vector[0], sizeof(double) * c->m2n0 * 2);
    dd_vector = c->d_vector[0];
  } else if (c->d_vector[1] == NULL) {
    cudaMalloc((void**)&c->d_vector[1], sizeof(double) * c->m2n0 * 2);
    dd_vector = c->d_vector[1];
  } else {
    std::swap(c->d_vector[0], c->d_vector[1]);
    dd_vector = c->d_vector[1];
  }
  cudaMemcpy(dd_vector, vector, sizeof(double) * c->m2n0 * 2, 
      cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] loading vector field data");

  /// 
  double *dd_jacobian;
  if (c->d_jacobian[0] == NULL) {
    cudaMalloc((void**)&c->d_jacobian[0], sizeof(double) * c->m2n0 * 4);
    dd_jacobian = c->d_jacobian[0];
  } else if (c->d_jacobian[1] == NULL) {
    cudaMalloc((void**)&c->d_jacobian[1], sizeof(double) * c->m2n0 * 4);
    dd_jacobian = c->d_jacobian[1];
  } else {
    std::swap(c->d_jacobian[0], c->d_jacobian[1]);
    dd_jacobian = c->d_jacobian[1];
  }
  cudaMemcpy(dd_jacobian, jacobian, sizeof(double) * c->m2n0 * 4, 
      cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] loading jacobian field data");
}

void xft_load_mesh(ctx_t *c,
    int nphi, int iphi, int vphi,
    int m2n0, int m2n1, int m2n2,
    const double *m2coords, const int *m2edges, const int *m2tris)
{
  c->nphi = nphi; 
  c->iphi = iphi;
  c->vphi = vphi;
  c->m2n0 = m2n0;
  c->m2n1 = m2n1;
  c->m2n2 = m2n2;

  cudaMalloc((void**)&c->d_m2coords, m2n0 * sizeof(double) * m2n0);
  cudaMemcpy(c->d_m2coords, m2coords, m2n0 * sizeof(double) * m2n0, cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_m2edges, m2n1 * sizeof(int) * m2n1);
  cudaMemcpy(c->d_m2edges, m2edges, m2n1 * sizeof(int) * m2n1, cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_m2tris, m2n2 * sizeof(int) * m2n2);
  cudaMemcpy(c->d_m2tris, m2tris, m2n2 * sizeof(int) * m2n2, cudaMemcpyHostToDevice);
  
  checkLastCudaError("[FTK-CUDA] loading xgc mesh");
}

void xft_load_interpolants(ctx_t *c, const std::vector<std::vector<xgc_interpolant_t>> &interpolants)
{
  assert(c->vphi == interpolants.size());

  cudaMalloc((void**)&c->d_interpolants, 
      c->m2n0 * sizeof(xgc_interpolant_t) * c->vphi);
  cudaMalloc((void**)&c->d_ptr_interpolants, c->vphi);

  std::vector<xgc_interpolant_t*> ptr_all_interpolants;

  for (int i = 0; i < interpolants.size(); i ++) {
    if (i == 0) ptr_all_interpolants.push_back(c->d_interpolants);
    else {
      ptr_all_interpolants.push_back(c->d_interpolants + (i-1) * c->m2n0);
      cudaMemcpy(c->d_interpolants + (i-1) * c->m2n0 * sizeof(xgc_interpolant_t), 
          interpolants[i].data(), c->m2n0 * sizeof(xgc_interpolant_t), cudaMemcpyHostToDevice);
    }
  }
  
  checkLastCudaError("[FTK-CUDA] loading xgc interpolants");
}
