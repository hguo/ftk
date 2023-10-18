#include <nvfunctional>
#include <ftk/numeric/mpas.hh>
#include <ftk/numeric/wachspress_interpolation.hh>
#include <ftk/filters/mpas_ocean_particle_tracker.cuh>
#include "common.cuh"

typedef mop_ctx_t ctx_t;

static const int MAX_VERTS = 10;
static const int MAX_LAYERS = 100;
static const double R0 = 6371229.0;

// what we need:
// - velocity field
// - vertical velocity field
// - zTop field
// - scalar attribute fields

__device__ __host__
inline bool point_in_cell(
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
  const int nverts = nedges_on_cell[cell];
  // double xv[MAX_VERTS][3];

  for (int i = 0; i < nverts; i ++) {
    for (int k = 0; k < 3; k ++) {
      const int vertex = verts_on_cell[cell * nverts + i];
      iv[k] = vertex;
      xv[i][k] = Xv[vertex*3+k];
    }
  }

  return ftk::point_in_mpas_cell<double>(nverts, xv, x);
}

__device__ __host__
static int locate_cell_local( // local search among neighbors
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
    for (int i = 0; i < nedges_on_cell[curr]; i ++) {
      const int cell = cells_on_cell[i + max_edges * curr];
      if (point_in_cell(
            cell, x, iv, xv,
            max_edges, Xv, nedges_on_cell, verts_on_cell))
        return cell;
    }
    return -1; // not found among neighbors
  }
}

__device__ __host__ 
static bool mpas_eval(
    const double *x,    // location
    double *v,          // return velocity
    double *vv,         // vertical velocity
    double *f,          // scalar attributs
    const double *V,    // velocity field
    const double *Vv,   // vertical velocities
    const double *zTop, // top layer depth
    const int nch,   // number of scalar attributes
    const double *A,    // scalar attributes
    const double *Xv,   // vertex locations
    const int max_edges,
    const int *nedges_on_cell, 
    const int *cells_on_cell,
    const int *verts_on_cell,
    const int nlayers,
    int &hint_c, 
    int &hint_l)        // hint for searching cell and layer
{
  int iv[MAX_VERTS];
  double xv[MAX_VERTS][3];

  const int cell = locate_cell_local(hint_c, 
      x, iv, xv, 
      Xv, max_edges, nedges_on_cell, 
      cells_on_cell, verts_on_cell);
  if (cell < 0) return false;
  else hint_c = cell;

  const int nverts = nedges_on_cell[cell];

  // compute weights based on xyzVerts
  double omega[MAX_VERTS]; 
  ftk::wachspress_weights(nverts, xv, x, omega); 

  // locate layer
  int layer = hint_l;
  double upper = 0.0, lower = 0.0;
  const double R = ftk::vector_2norm<3>(x);
  const double z = R - R0;
  int dir; // 0=up, 1=down
 
  bool succ = false;
  if (layer >= 0) { // interpolate upper/lower tops and check if x remains in the layer
    layer = hint_l;
    for (int i = 0; i < nverts; i ++) {
      upper += omega[i] * zTop[ iv[i] * nlayers + layer ];
      lower += omega[i] * zTop[ iv[i] * nlayers + layer+1 ];

      if (z > upper)
        dir = 0; // up
      else if (z <= lower) 
        dir = 1; // down
      else 
        succ = true;
    }
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

    for (int k = 0; k < nch; k ++)
      f[k] += omega[i] * (
                alpha * A[ k + nch * (iv[i] * nlayers + layer) ]
              + beta  * A[ k + nch * (iv[i] * nlayers + layer + 1) ]);

    *vv +=   alpha * Vv[ iv[i] * (nlayers + 1) + layer ]
           + beta  * Vv[ iv[i] * (nlayers + 1) + layer + 1];
  }

  return true;
}

///////////////////////////
void mop_create_ctx(mop_ctx_t **c_, int device)
{
  *c_ = (ctx_t*)malloc(sizeof(ctx_t));
  ctx_t *c = *c_;
  memset(c, 0, sizeof(ctx_t));

  c->device = device;
  cudaSetDevice(device);

  c->d_Xc = NULL;
  c->d_Xv = NULL;
  c->d_nedges_on_cell = NULL;
  c->d_cells_on_cell = NULL;
  c->d_verts_on_cell = NULL;
 
  c->d_V[0] = NULL;
  c->d_V[1] = NULL;
}

void mop_destroy_ctx(mop_ctx_t **c_)
{
  ctx_t *c = *c_;

  if (c->d_Xc != NULL) cudaFree(c->d_Xc);
  if (c->d_Xv != NULL) cudaFree(c->d_Xv);
  if (c->d_nedges_on_cell != NULL) cudaFree(c->d_nedges_on_cell);
  if (c->d_cells_on_cell != NULL) cudaFree(c->d_cells_on_cell);
  if (c->d_verts_on_cell != NULL) cudaFree(c->d_verts_on_cell);

  if (c->d_V[0] != NULL) cudaFree(c->d_V[0]);
  if (c->d_V[1] != NULL) cudaFree(c->d_V[1]);

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
    const int nch,
    const double *Xc,
    const double *Xv,
    const int *nedges_on_cell, 
    const int *cells_on_cell,
    const int *verts_on_cell)
{
  c->ncells = ncells;
  c->nlayers = nlayers;
  c->nverts = nverts;
  c->max_edges = max_edges;
  c->nch = nch;

  cudaMalloc((void**)&c->d_Xc, size_t(ncells) * sizeof(double) * 3);
  cudaMemcpy(c->d_Xc, Xc, size_t(ncells) * sizeof(double) * 3, cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_Xv, size_t(nverts) * sizeof(double) * 3);
  cudaMemcpy(c->d_Xv, Xv, size_t(nverts) * sizeof(double) * 3, cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_nedges_on_cell, size_t(ncells) * sizeof(int));
  cudaMemcpy(c->d_nedges_on_cell, nedges_on_cell, size_t(ncells) * sizeof(int), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_cells_on_cell, size_t(ncells) * max_edges * sizeof(int));
  cudaMemcpy(c->d_cells_on_cell, cells_on_cell, size_t(ncells) * max_edges * sizeof(int), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_verts_on_cell, size_t(nverts) * 3 * sizeof(int));
  cudaMemcpy(c->d_verts_on_cell, verts_on_cell, size_t(nverts) * 3 * sizeof(int), cudaMemcpyHostToDevice);

  checkLastCudaError("[FTK-CUDA] loading mpas mesh");
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

void mop_load_data(mop_ctx_t *c, const double *V) 
{
  fprintf(stderr, "nch=%d, ncells=%d, nlayers=%d\n", c->nch, c->ncells, c->nlayers);
  load_data<double>(c->d_V, V, size_t(c->nch) * c->ncells * c->nlayers, "V");
}


void mop_execute(mop_ctx_t *c, int current_timestep)
{
  size_t ntasks = c->nparticles;
  fprintf(stderr, "ntasks=%zu\n", ntasks);
  
  const int maxGridDim = 1024;
  const int blockSize = 256;
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;

  if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else gridSize = dim3(nBlocks);

  // TODO
#if 0
  sweep_simplices<int, double><<<gridSize, blockSize>>>(
      scope, current_timestep, 
      c->factor,
      c->nphi, c->iphi, c->vphi, 
      c->m2n0, c->m2n1, c->m2n2, 
      c->d_m2coords, c->d_m2edges, c->d_m2tris, 
      c->d_psin,
      c->d_interpolants, 
      c->d_V[0], c->d_V[1],
      c->d_vector[0], c->d_vector[1],
      c->d_jacobian[0], c->d_jacobian[1], 
      *c->dncps, c->dcps);
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] sweep_simplicies");

  cudaMemcpy(&c->hncps, c->dncps, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
  cudaMemset(c->dncps, 0, sizeof(unsigned long long)); // clear the counter
  checkLastCudaError("[FTK-CUDA] cuda memcpy device to host, 1");
  fprintf(stderr, "ncps=%llu\n", c->hncps);
  cudaMemcpy(c->hcps, c->dcps, sizeof(cp_t) * c->hncps, cudaMemcpyDeviceToHost);
  
  checkLastCudaError("[FTK-CUDA] cuda memcpy device to host, 2");
#endif
}

void mop_swap(mop_ctx_t *c)
{
  std::swap(c->d_V[0], c->d_V[1]);
}
