#include <nvfunctional>
#include <ftk/numeric/mpas.hh>
#include <ftk/numeric/wachspress_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/filters/mpas_ocean_particle_tracker.cuh>
#include <ftk/numeric/rk4.hh>
#include "common.cuh"

typedef mop_ctx_t ctx_t;

static const int MAX_VERTS = 10;
static const int MAX_ATTRS = 10;
// static const int MAX_LAYERS = 100;
static const double R0 = 6371229.0;
  
static const int maxGridDim = 1024;
static const int blockSize = 256;

// what we need:
// - velocity field
// - vertical velocity field
// - zTop field
// - scalar attribute fields

__device__ __host__
inline static int c2ci(int c) { return c - 1; }

__device__ __host__
inline static int v2vi(int v) { return v - 1; }

__device__ __host__
inline static bool point_in_cell(
    const int cell,
    const double *x,
    int iv[],
    double xv[][3], // returns vertex coordinates
    const int max_edges,
    const double *Xv,
    const int *nedges_on_cell, 
    const int *verts_on_cell)
{
  // if (cell < 0) return false;
  const int nverts = nedges_on_cell[c2ci(cell)];
  // double xv[MAX_VERTS][3];

  for (int i = 0; i < nverts; i ++) {
    const int vertex = verts_on_cell[c2ci(cell) * max_edges + i];
    iv[i] = v2vi(vertex);
    for (int k = 0; k < 3; k ++)
      xv[i][k] = Xv[iv[i]*3+k];
  }

  bool succ = ftk::point_in_mpas_cell<double>(nverts, xv, x);
#if 0
  printf("x=%f, %f, %f, cell=%d, nverts=%d, max_edges=%d, iv=%d, %d, %d, %d, %d, %d, %d, succ=%d\n", 
      x[0], x[1], x[2], cell, nverts, max_edges,
      iv[0], iv[1], iv[2], iv[3], iv[4], iv[5], iv[6], succ);
#endif
  return succ;
}

__device__ __host__
inline static int locate_cell_local( // local search among neighbors
    const int curr, // current cell
    const double *x,
    int iv[], // returns vertex ids
    double xv[][3], // returns vertex coordinates
    const double *Xv,
    const int max_edges,
    const int *nedges_on_cell, // also n_verts_on_cell
    const int *cells_on_cell,
    const int *verts_on_cell)
{
  if (curr < 0)
    return -1; // not found
  else if (point_in_cell(
        curr, x, iv, xv, 
        max_edges, Xv, nedges_on_cell, verts_on_cell))
    return curr;
  else {
    for (int i = 0; i < nedges_on_cell[c2ci(curr)]; i ++) {
      const int cell = cells_on_cell[i + max_edges * c2ci(curr)];
      if (point_in_cell(
            cell, x, iv, xv,
            max_edges, Xv, nedges_on_cell, verts_on_cell)) {
        // printf("curr=%d, cell=%d\n", curr, cell);
        return cell;
      }
    }
    return -1; // not found among neighbors
  }
}

__device__ __host__ 
inline static bool mpas_eval_static(
    const double *x,    // location
    double *v,          // return velocity
    double *vv,         // vertical velocity
    double *f,          // scalar attributs
    const double *V,    // velocity field
    const double *Vv,   // vertical velocities
    const double *zTop, // top layer depth
    const int attrs,    // number of scalar attributes
    const double *A,    // scalar attributes
    const double *Xv,   // vertex locations
    const int max_edges,
    const int *nedges_on_cell, 
    const int *cells_on_cell,
    const int *verts_on_cell,
    const int nlayers,
    unsigned long long &hint_c, 
    unsigned int &hint_l)        // hint for searching cell and layer
{
  int iv[MAX_VERTS];
  double xv[MAX_VERTS][3];

#if 0
  printf("hc=%llu, hl=%lu, x=%f, %f, %f\n", 
      hint_c, hint_l, x[0], x[1], x[2]);
  printf("c=%llu, l=%lu, x=%f, %f, %f, v=%f, %f, %f, %f\n", 
      hint_c, hint_l, 
      x[0], x[1], x[2], 
      v[0], v[1], v[2], v[3]);
#endif

  const int cell = locate_cell_local(hint_c, 
      x, iv, xv, 
      Xv, max_edges, nedges_on_cell, 
      cells_on_cell, verts_on_cell);
  if (cell < 0) return false;
  else hint_c = cell;

  const int nverts = nedges_on_cell[c2ci(cell)];

  // compute weights based on xyzVerts
  double omega[MAX_VERTS]; 
  ftk::wachspress_weights(nverts, xv, x, omega); 

#if 0
  printf("cell=%d, nverts=%d, x=%f, %f, %f, iv=%d, %d, %d, %d, %d, %d, %d, weights=%f, %f, %f, %f, %f, %f, %f\n", 
      cell, nverts,
      x[0], x[1], x[2],
      iv[0], iv[1], iv[2], iv[3], iv[4], iv[5], iv[6],
      omega[0], omega[1], omega[2], omega[3], omega[4], omega[5], omega[6]);
#endif

  // locate layer
  int layer = hint_l;
  double upper = 0.0, lower = 0.0;
  const double R = ftk::vector_2norm<3>(x);
  const double z = R - R0;
  int dir; // 0=up, 1=down

  // printf("hint_l=%d\n", hint_l);

  bool succ = false;
  if (layer >= 0) { // interpolate upper/lower tops and check if x remains in the layer
    layer = hint_l;
    for (int i = 0; i < nverts; i ++) {
      // printf("%d, %d, %d\n", iv[i], nlayers, layer);
      upper += omega[i] * zTop[ iv[i] * nlayers + layer ];
      lower += omega[i] * zTop[ iv[i] * nlayers + layer+1 ];
    }

    // printf("z=%f, lower=%f, upper=%f\n", z, lower, upper);

    if (z > upper)
      dir = 0; // up
    else if (z <= lower) 
      dir = 1; // down
    else 
      succ = true;

  } else {
    layer = 0;
    dir = 1;
  }

  if (!succ) {
    if (dir == 1) { // downward
      upper = lower;
      for (layer = layer + 1 ; layer < nlayers-1; layer ++) {
        lower = 0.0;
        for (int k = 0; k < nverts; k ++)
          lower += omega[k] * zTop[ iv[k] * nlayers + layer + 1];

        if (z <= upper && z > lower) {
          succ = true;
          break;
        } else 
          upper = lower;
      }
    } else { // upward
      lower = upper;
      for (layer = layer - 1; layer >= 0; layer --) {
        upper = 0.0;
        for (int k = 0; k < nverts; k ++)
          upper += omega[k] * zTop[ iv[k] * nlayers + layer];

        if (z <= upper && z > lower) {
          succ = true;
          break;
        } else 
          lower = upper;
      }
    }
  }

  if (!succ) 
    return false;

  hint_l = layer;

  const double alpha = (z - lower) / (upper - lower), 
               beta = 1.0 - alpha;

  // reset values before interpolation
  memset(v, 0, sizeof(double)*3);
  memset(f, 0, sizeof(double)*3);
  *vv = 0.0;

  // interpolation
  for (int i = 0; i < nverts; i ++) {
    for (int k = 0; k < 3; k ++)
      v[k] += omega[i] * (
                alpha * V[ k + 3 * (iv[i] * nlayers + layer) ]
              + beta  * V[ k + 3 * (iv[i] * nlayers + layer + 1) ]);

    for (int k = 0; k < attrs; k ++)
      f[k] += omega[i] * (
                alpha * A[ k + attrs * (iv[i] * nlayers + layer) ]
              + beta  * A[ k + attrs * (iv[i] * nlayers + layer + 1) ]);

    *vv +=   alpha * Vv[ iv[i] * (nlayers + 1) + layer ]
           + beta  * Vv[ iv[i] * (nlayers + 1) + layer + 1];
  }

  return true;
}

__device__ __host__ 
inline static bool mpas_eval(
    const double *x,             // location
    double *v,                   // return velocity
    double *vv,                  // vertical velocity
    double *f,                   // scalar attributs
    const double *const V[2],    // velocity field
    const double *const Vv[2],   // vertical velocities
    const double *const zTop[2], // top layer depth
    const int attrs,             // number of scalar attributes
    const double *const A[2],    // scalar attributes
    const double *Xv,            // vertex locations
    const int max_edges,
    const int *nedges_on_cell, 
    const int *cells_on_cell,
    const int *verts_on_cell,
    const int nlayers,
    unsigned long long &hint_c, 
    unsigned int &hint_l)        // hint for searching cell and layer
{
  if (V[1]) { // two timesteps
    double _v[2][3], _vv[2], _f[2][MAX_ATTRS];

    for (int i = 0; i < 2; i ++) 
      if (!mpas_eval_static(x, 
          _v[i], &_vv[i], _f[i], 
          V[i], Vv[i], zTop[i], attrs, A[i], 
          Xv, max_edges, nedges_on_cell, cells_on_cell, verts_on_cell, 
          nlayers, hint_c, hint_l))
        return false;

  } else
    return mpas_eval_static(
        x, v, vv, f, V[0], Vv[0], zTop[0], attrs, A[0], 
        Xv, max_edges, nedges_on_cell, cells_on_cell, verts_on_cell, 
        nlayers, hint_c, hint_l);
}

template <int ORDER=1>
__device__
inline static bool spherical_rk_with_vertical_velocity(
    const double h,
    ftk::feature_point_lite_t& p, 
    double *v0,                  // return velocity
    double *vv,                  // vertical velocity
    double *f,                   // scalar attributs
    const double *const V[2],    // velocity field
    const double *const Vv[2],   // vertical velocities
    const double *const zTop[2], // top layer depth
    const int nattrs,            // number of scalar attributes
    const double *const A[2],    // scalar attributes
    const double *const Xv,      // vertex locations
    const int max_edges,
    const int *nedges_on_cell, 
    const int *cells_on_cell,
    const int *verts_on_cell,
    const int nlayers,
    unsigned long long &hint_c, 
    unsigned int &hint_l)        // hint for searching cell and layer
{

  if (ORDER == 1) {
    double v[4];
    if (!mpas_eval(p.x, v, vv, f,
          V, Vv, zTop, nattrs, A, Xv, 
          max_edges, nedges_on_cell, cells_on_cell, verts_on_cell, 
          nlayers, hint_c, hint_l))
      return false;

    for (int k = 0; k < 4; k ++)
      v0[k] = v[k];

    ftk::angular_stepping(p.x, v, h, p.x);
    return true;
  } else if (ORDER == 4) {
    return true;
  } else
    return false;
  
  // printf("x=%f, %f, %f, v=%f, %f, %f, %f\n", 
  //     p.x[0], p.x[1], p.x[2], v[0], v[1], v[2], v[3]);
}

__global__
static void initialize_c2v(
    const int ncells,
    const int nverts,
    const double *Xc, // cell x
    const double *Xv, // vertex x
    const int *cells_on_vert,
    double *interpolants,
    bool *cell_on_boundary)
{
  unsigned long long i = getGlobalIdx_3D_1D();
  if (i >= nverts) return;

  bool boundary = false;
  double xc[3][3], xv[3];

  const int vi = i;
#if 0
  printf("%d, %d, %d, %d, %d, %d, %d\n",
      vi, 
      cells_on_vert[vi*3],
      cells_on_vert[vi*3+1],
      cells_on_vert[vi*3+2], 
      vi*3, vi*3+1, vi*3+2);
#endif

  for (int j = 0; j < 3; j ++) {
    const int ci = c2ci( cells_on_vert[j + vi*3] );
    // const int ci = cells_on_vert[j + i*3];
    // printf("%d, %d\n", vi, ci);
    if (ci < 0) {
      boundary = true;
      break;
    }
    xv[j] = Xv[j + vi*3];
    for (int k = 0; k < 3; k ++)
      xc[j][k] = Xc[k + ci*3];
  }

  if (boundary) {
    cell_on_boundary[vi] = true;
  } else {
    cell_on_boundary[vi] = false;
    
    double mu[3];
    bool succ = ftk::inverse_lerp_s2v3_0(xc, xv, mu);
    for (int k = 0; k < 3; k ++)
      interpolants[k + vi*3] = mu[k];
  }
}


__global__
static void interpolate_c2v(
    double *V, // vertexwise data
    const double *C, // cellwise data
    const int nch,
    const int nlayers,
    const int ncells,
    const int nverts,
    const double *interpolants,
    const int *cells_on_vert,
    const bool *cell_on_boundary)
{
  // unsigned long long i = 192;
  // if (threadIdx.x > 0) return;
  unsigned long long i = getGlobalIdx_3D_1D();
  if (i >= nverts) return;

  const int vi = i;

  const bool boundary = cell_on_boundary[vi];
  double lambda[3];
  int ci[3];

  for (int l = 0; l < 3; l ++) {
    lambda[l] = interpolants[l + vi * 3];
    ci[l] = c2ci(cells_on_vert[l + vi * 3]);
  }

  // printf("ci=%d, %d, %d, boundary=%d, lambda=%f, %f, %f\n", 
  //     ci[0], ci[1], ci[2], boundary, 
  //     lambda[0], lambda[1], lambda[2]);

  for (int layer = 0; layer < nlayers; layer ++) {
    for (int k = 0; k < nch; k ++) {
      double &v = V[k + nch * (layer + vi * nlayers)];
      if (boundary) {
        v = 0;
      } else {
        bool invalid = false;
        for (int l = 0; l < 3; l ++) {
          double val = C[k + nch * (layer + ci[l] * nlayers)];
          if (val < -1e33) {
            invalid = true;
            break;
          }
          v += lambda[l] * val;
        }
        
        if (invalid)
          v = nan("1");
      }
    }
  }
}

__global__
static void mpas_trace(
    const int nsteps,
    const double h,
    const int nparticles,
    ftk::feature_point_lite_t* particles,
    const double *const V[2],    // velocity field
    const double *const Vv[2],   // vertical velocities
    const double *const zTop[2], // top layer depth
    const int nattrs,   // number of scalar attributes
    const double *const A[2],    // scalar attributes
    const double *Xv,   // vertex locations
    const int max_edges,
    const int *nedges_on_cell, 
    const int *cells_on_cell,
    const int *verts_on_cell,
    const int nlayers)
{
  unsigned long long i = getGlobalIdx_3D_1D();
  if (i >= nparticles) return;

  ftk::feature_point_lite_t &p = particles[i];
  double v0[3], vv[3], f[MAX_ATTRS];

  unsigned long long &hint_c = p.tag;
  unsigned int &hint_l = p.type;

  for (int j = 0; j < nsteps; j ++) {
    bool succ = spherical_rk_with_vertical_velocity<1>(
        h, p, v0, vv, f, 
        V, Vv, zTop, nattrs, A, Xv, 
        max_edges, nedges_on_cell, cells_on_cell, verts_on_cell, 
        nlayers, hint_c, hint_l);

    if (!succ) 
      break;
  }
}

///////////////////////////
void mop_create_ctx(mop_ctx_t **c_, int device)
{
  *c_ = (ctx_t*)malloc(sizeof(ctx_t));
  ctx_t *c = *c_;
  memset(c, 0, sizeof(ctx_t));

  c->device = device;
  cudaSetDevice(device);
}

void mop_destroy_ctx(mop_ctx_t **c_)
{
  ctx_t *c = *c_;
  if (c == NULL) return;

  if (c->d_Xc != NULL) cudaFree(c->d_Xc);
  if (c->d_Xv != NULL) cudaFree(c->d_Xv);
  if (c->d_nedges_on_cell != NULL) cudaFree(c->d_nedges_on_cell);
  if (c->d_cells_on_cell != NULL) cudaFree(c->d_cells_on_cell);
  if (c->d_verts_on_cell != NULL) cudaFree(c->d_verts_on_cell);
  if (c->d_cells_on_vert != NULL) cudaFree(c->d_cells_on_vert);

  for (int i = 0; i < 2; i ++) {
    if (c->d_V[i] != NULL) cudaFree(c->d_V[i]);
    if (c->d_Vv[i] != NULL) cudaFree(c->d_Vv[i]);
    if (c->d_zTop[i] != NULL) cudaFree(c->d_zTop[i]);
    if (c->d_A[i] != NULL) cudaFree(c->d_A[i]);
  }

  if (c->dd_V) cudaFree(c->dd_V);
  if (c->dd_Vv) cudaFree(c->dd_Vv);
  if (c->dd_zTop) cudaFree(c->dd_zTop);
  if (c->dd_A) cudaFree(c->dd_A);

  if (c->d_c2v_interpolants != NULL) cudaFree(c->d_c2v_interpolants);
  if (c->d_cell_on_boundary != NULL) cudaFree(c->d_cell_on_boundary);
  if (c->dcw != NULL) cudaFree(c->dcw);

  if (c->dparts != NULL) cudaFree(c->dparts);
  if (c->hparts != NULL) free(c->hparts);

  free(*c_);
  *c_ = NULL;
}

template <typename T=double>
static void realloc_both(
    T** hbuf,
    T** dbuf,
    const size_t n, 
    const size_t m)
{
  if (*hbuf == NULL)
    *hbuf = (T*)malloc(m * sizeof(T));
  else if (m != n)
    *hbuf = (T*)realloc(hbuf, m * sizeof(T));

  if (*dbuf == NULL)
    cudaMalloc((void**)dbuf, m * sizeof(T));
  else if (m != n) {
    cudaFree(dbuf);
    cudaMalloc((void**)dbuf, m * sizeof(T));
  }
}

void mop_load_particles(mop_ctx_t *c, 
    const int n, 
    ftk::feature_point_lite_t *buf)
{
  realloc_both(&c->hparts, &c->dparts, c->nparticles, n);
  c->nparticles = n;

  cudaMemcpy(c->dparts, buf, n * sizeof(ftk::feature_point_lite_t), 
      cudaMemcpyHostToDevice);
  
  checkLastCudaError("[FTK-CUDA] load particles");
}

void mop_load_mesh(mop_ctx_t *c,
    const int ncells, 
    const int nlayers, 
    const int nverts, 
    const int max_edges,
    const int nattrs,
    const double *Xc,
    const double *Xv,
    const int *nedges_on_cell, 
    const int *cells_on_cell,
    const int *verts_on_cell,
    const int *cells_on_vert)
{
  c->ncells = ncells;
  c->nlayers = nlayers;
  c->nverts = nverts;
  c->max_edges = max_edges;
  c->nattrs = nattrs;

  cudaMalloc((void**)&c->d_Xc, 
      size_t(ncells) * sizeof(double) * 3);
  cudaMemcpy(c->d_Xc, Xc, 
      size_t(ncells) * sizeof(double) * 3, 
      cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_Xv, 
      size_t(nverts) * sizeof(double) * 3);
  cudaMemcpy(c->d_Xv, Xv, 
      size_t(nverts) * sizeof(double) * 3, 
      cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_nedges_on_cell, 
      size_t(ncells) * sizeof(int));
  cudaMemcpy(c->d_nedges_on_cell, nedges_on_cell, 
      size_t(ncells) * sizeof(int), 
      cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_cells_on_cell, 
      size_t(ncells) * max_edges * sizeof(int));
  cudaMemcpy(c->d_cells_on_cell, cells_on_cell, 
      size_t(ncells) * max_edges * sizeof(int), 
      cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_verts_on_cell, 
      size_t(nverts) * max_edges * sizeof(int));
  cudaMemcpy(c->d_verts_on_cell, verts_on_cell, 
      size_t(nverts) * max_edges * sizeof(int), 
      cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_cells_on_vert, 
      size_t(nverts) * 3 * sizeof(int));
  cudaMemcpy(c->d_cells_on_vert, cells_on_vert, 
      size_t(nverts) * 3 * sizeof(int),
      cudaMemcpyHostToDevice);

  checkLastCudaError("[FTK-CUDA] loading mpas mesh");

  // initialize interpolants
  {
    cudaMalloc((void**)&c->d_c2v_interpolants, 
        size_t(c->nverts) * 3 * sizeof(double));

    cudaMalloc((void**)&c->d_cell_on_boundary,
        size_t(c->ncells) * sizeof(bool));
    
    // cudaDeviceSynchronize();
    checkLastCudaError("[FTK-CUDA] initializing mpas c2v interpolants: malloc");

    size_t ntasks = c->nverts;
    const int nBlocks = idivup(ntasks, blockSize);
    dim3 gridSize;
  
    if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
    else gridSize = dim3(nBlocks);

    initialize_c2v<<<gridSize, blockSize>>>(
        c->ncells, 
        c->nverts, 
        c->d_Xc,
        c->d_Xv,
        c->d_cells_on_vert, 
        c->d_c2v_interpolants,
        c->d_cell_on_boundary);
 
    fprintf(stderr, "c2v initialization done.\n");
    // cudaDeviceSynchronize();
    checkLastCudaError("[FTK-CUDA] initializing mpas c2v interpolants");
  }
}

template <typename T=double>
static void load_data(
    T** dbuf,
    const T *buf,
    const size_t n, 
    const char *name)
{
  T *d;
  if (dbuf[0] == NULL) {
    cudaMalloc((void**)&dbuf[0], sizeof(T) * n);
    checkLastCudaError("[FTK-CUDA] loading data 0");
    d = dbuf[0];
  } else if (dbuf[1] == NULL) {
    cudaMalloc((void**)&dbuf[1], sizeof(T) * n);
    checkLastCudaError("[FTK-CUDA] loading data 1");
    d = dbuf[1];
  } else {
    std::swap(dbuf[0], dbuf[1]);
    d = dbuf[1];
  }
  cudaMemcpy(d, buf, sizeof(T) * n, cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] cudaMemcpy");
}

static void load_cw_data(
    mop_ctx_t *c,
    double **dbuf, // an array with two device pointers
    double ***ddbuf, // device array of pointers
    const double *cw,
    const int nch,
    const int nlayers) // host pointer with cw data
{
  // initialize buffers for vertexwise data in two adjacent timesteps
  double *d;
  if (dbuf[0] == NULL) {
    cudaMalloc((void**)&dbuf[0], 
        sizeof(double) * c->nverts * nch * nlayers);
    d = dbuf[0];
    checkLastCudaError("malloc vw buffers 0");
  } else if (dbuf[1] == NULL) {
    cudaMalloc((void**)&dbuf[1], 
        sizeof(double) * c->nverts * nch * nlayers);
    d = dbuf[1];
    checkLastCudaError("malloc vw buffers 1");
  } else {
    std::swap(dbuf[0], dbuf[1]);
    d = dbuf[1];
  }
 
  // cudaDeviceSynchronize();
  // checkLastCudaError("memcpy to c2w buffer: dev ptrs0");
  if (*ddbuf == NULL) {
    cudaMalloc((void**)ddbuf, sizeof(double*) * 2);
    checkLastCudaError("memcpy to c2w buffer: allocating ddbuf");
  }
  // fprintf(stderr, "%p\n", ddbuf);
  cudaMemcpy(*ddbuf, dbuf, sizeof(double*) * 2, 
      cudaMemcpyHostToDevice);
  checkLastCudaError("memcpy to c2w buffer: dev ptrs");
  
  // copy data to c2w buffer
  assert(c->dcw != NULL);
  cudaMemcpy(c->dcw, cw, 
      sizeof(double) * c->ncells * nch * nlayers, 
      cudaMemcpyHostToDevice);
  // cudaDeviceSynchronize();
  // fprintf(stderr, "dcw=%p\n", c->dcw);
  checkLastCudaError("memcpy to c2w buffer");

  // c2w interpolation
  {
    size_t ntasks = c->nverts;
    const int nBlocks = idivup(ntasks, blockSize);
    dim3 gridSize;
  
    if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
    else gridSize = dim3(nBlocks);

    interpolate_c2v<<<gridSize, blockSize>>>(
        d, c->dcw, nch, nlayers,
        c->ncells,
        c->nverts,
        c->d_c2v_interpolants,
        c->d_cells_on_vert,
        c->d_cell_on_boundary);
 
    // cudaDeviceSynchronize();
    checkLastCudaError("c2w interpolation");
  }
}

void mop_load_data(mop_ctx_t *c, 
    const double *V,
    const double *Vv,
    const double *zTop,
    const double *A)
{
  load_data<double>(c->d_V, V, 3 * c->nverts * c->nlayers, "V");
  load_data<double>(c->d_Vv, Vv, c->nverts * (c->nlayers+1), "Vv");
  load_data<double>(c->d_zTop, zTop, c->nverts * c->nlayers, "zTop");
  load_data<double>(c->d_A, A, c->nattrs * c->nverts * c->nlayers, "f");
}

void mop_load_data_cw(mop_ctx_t *c,
    const double *Vc,
    const double *Vvc,
    const double *zTopc,
    const double *Ac)
{
  if (c->dcw == NULL)
    cudaMalloc((void**)&c->dcw, 
        sizeof(double) * c->ncells * 3 * (c->nlayers+1)); // a sufficiently large buffer for loading cellwise data
  // cudaDeviceSynchronize();
  checkLastCudaError("malloc dcw");

  load_cw_data(c, c->d_V, &c->dd_V, Vc, 3, c->nlayers); // V
  load_cw_data(c, c->d_Vv, &c->dd_Vv, Vvc, 1, c->nlayers+1); // Vv
  load_cw_data(c, c->d_zTop, &c->dd_zTop, zTopc, 1, c->nlayers); // zTop
  load_cw_data(c, c->d_A, &c->dd_A, Ac, c->nattrs, c->nlayers); // f
}

void mop_execute(mop_ctx_t *c, int current_timestep)
{
  size_t ntasks = c->nparticles;
  // fprintf(stderr, "ntasks=%zu\n", ntasks);
  
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;

  if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else gridSize = dim3(nBlocks);

  // fprintf(stderr, "gridSize=%d, %d, %d, blockSize=%d\n", gridSize.x, gridSize.y, gridSize.z, blockSize);
  // cudaDeviceSynchronize();
  // checkLastCudaError("[FTK-CUDA] mop_execute, pre");

  mpas_trace<<<gridSize, blockSize>>>(
      32768, 0.0001, 
      ntasks, 
      c->dparts,
      c->dd_V,
      c->dd_Vv,
      c->dd_zTop, 
      c->nattrs, 
      c->dd_A,
      c->d_Xv,
      c->max_edges, 
      c->d_nedges_on_cell, 
      c->d_cells_on_cell,
      c->d_verts_on_cell, 
      c->nlayers);
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] mop_execute");

  cudaMemcpy(c->hparts, c->dparts, 
      c->nparticles * sizeof(ftk::feature_point_lite_t), 
      cudaMemcpyDeviceToHost);
  checkLastCudaError("[FTK-CUDA] mop_execute: memcpy");

  fprintf(stderr, "exiting kernel\n");
}

void mop_swap(mop_ctx_t *c)
{
  std::swap(c->d_V[0], c->d_V[1]);
}
