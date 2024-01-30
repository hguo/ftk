#include <nvfunctional>
#include <ftk/numeric/mpas.hh>
#include <ftk/numeric/wachspress_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/filters/mpas_ocean_particle_tracker.cuh>
#include <ftk/numeric/rk4.hh>
#include <ftk/numeric/rad.hh>
#include <ftk/basic/kd_lite.hh>
#include <assert.h>
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

// template type
// - F, float type for computing
// - Fm, float type for mesh coordinates
// - Fv, float type for variables

__device__ __host__
inline static int c2ci(int c) { return c - 1; }

__device__ __host__
inline static int e2ei(int e) { return e - 1; }

__device__ __host__
inline static int v2vi(int v) { return v - 1; }

template <typename F=double, typename Fm=double>
__device__ __host__
inline static bool point_in_cell(
    const int cell,
    const F *x,
    int iv[],
    F xv[][3], // returns vertex coordinates
    const int max_edges,
    const Fm *Xv,
    const int *nedges_on_cell, 
    const int *verts_on_cell)
{
  // if (cell < 0) return false;
  const int nverts = nedges_on_cell[c2ci(cell)];
  // F xv[MAX_VERTS][3];

  for (int i = 0; i < nverts; i ++) {
    const int vertex = verts_on_cell[c2ci(cell) * max_edges + i];
    iv[i] = v2vi(vertex);
    for (int k = 0; k < 3; k ++)
      xv[i][k] = Xv[iv[i]*3+k];
  }

  bool succ = ftk::point_in_mpas_cell<F>(nverts, xv, x);
#if 0
  printf("x=%f, %f, %f, cell=%d, nverts=%d, max_edges=%d, iv=%d, %d, %d, %d, %d, %d, %d, succ=%d\n", 
      x[0], x[1], x[2], cell, nverts, max_edges,
      iv[0], iv[1], iv[2], iv[3], iv[4], iv[5], iv[6], succ);
#endif
  return succ;
}

template <typename F=double, typename Fm=double>
__device__ __host__
static int locate_cell_global( // global cell search
    const F *x,
    int iv[], // returns vertex ids
    F xv[][3], // returns vertex coordinates
    const int ncells,
    const int *kdheap,
    const Fm *Xc,
    const Fm *Xv,
    const int max_edges,
    const int *nedges_on_cell, // also n_verts_on_cell
    const int *cells_on_cell,
    const int *verts_on_cell)
{
  const int c = ftk::kdlite_nearest<3, int, Fm>(
      ncells, Xc, kdheap, x) + 1; // to cell id
 
  if (point_in_cell<F, Fm>(
        c, x, iv, xv, 
        max_edges, Xv, nedges_on_cell, verts_on_cell)) {

    // printf("cell=%d, x=%f, %f, %f\n", c, x[0], x[1], x[2]);
    return c;
  }
  else 
    return -1;
}

template <typename F=double, typename Fm=double>
__device__ __host__
static int locate_cell_local( // local search among neighbors
    const int curr, // current cell
    const F *x,
    int iv[], // returns vertex ids
    F xv[][3], // returns vertex coordinates
    const Fm *Xv,
    const int max_edges,
    const int *nedges_on_cell, // also n_verts_on_cell
    const int *cells_on_cell,
    const int *verts_on_cell)
{
  if (curr < 0)
    return -1; // not found
  else if (point_in_cell<F, Fm>(
        curr, x, iv, xv, 
        max_edges, Xv, nedges_on_cell, verts_on_cell))
    return curr;
  else {
    for (int i = 0; i < nedges_on_cell[c2ci(curr)]; i ++) {
      const int cell = cells_on_cell[i + max_edges * c2ci(curr)];
      if (point_in_cell<F, Fm>(
            cell, x, iv, xv,
            max_edges, Xv, nedges_on_cell, verts_on_cell)) {
        // printf("curr=%d, cell=%d\n", curr, cell);
        return cell;
      }
    }
    return -1; // not found among neighbors
  }
}

template <typename F=double, typename Fm=double, typename Fv=double>
__device__ __host__ 
static bool mpas_eval_static(
    const F *x,    // location
    F *v,          // return velocity
    F *vv,         // vertical velocity
    F *f,          // scalar attributs
    const Fv *V,    // velocity field
    const Fv *Vv,   // vertical velocities
    const Fv *zTop, // top layer depth
    const int nattrs,    // number of scalar attributes
    const Fv *A,    // scalar attributes
    const Fm *Xv,   // vertex locations
    const int max_edges,
    const int *nedges_on_cell, 
    const int *cells_on_cell,
    const int *verts_on_cell,
    const int nlayers,
    unsigned long long &hint_c, 
    unsigned int &hint_l)        // hint for searching cell and layer
{
  int iv[MAX_VERTS];
  F xv[MAX_VERTS][3];

#if 0
  printf("hc=%llu, hl=%u, x=%f, %f, %f\n", 
      hint_c, hint_l, x[0], x[1], x[2]);
  printf("c=%llu, l=%u, x=%f, %f, %f, v=%f, %f, %f, %f\n", 
      hint_c, hint_l, 
      x[0], x[1], x[2], 
      v[0], v[1], v[2], v[3]);
#endif

  const int cell = locate_cell_local<F, Fm>(hint_c, 
      x, iv, xv, 
      Xv, max_edges, nedges_on_cell, 
      cells_on_cell, verts_on_cell);
  if (cell < 0) return false;
  else hint_c = cell;

  const int nverts = nedges_on_cell[c2ci(cell)];

  // compute weights based on xyzVerts
  F omega[MAX_VERTS]; 
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
  if (layer >= nlayers || layer < 0)
    layer = 0;

  F upper = 0.0, lower = 0.0;
  const F R = ftk::vector_2norm<3>(x);
  const F z = R - R0;
  int dir; // 0=up, 1=down

  // printf("hint_l=%d\n", hint_l);

  bool succ = false;
  if (layer >= 0) { // interpolate upper/lower tops and check if x remains in the layer
    // layer = hint_l;
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
        
        // printf("moving downward, layer=%d, z=%f, upper=%f, lower=%f\n", 
        //     layer, z, upper, lower);

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
        
        // printf("moving upward, layer=%d, z=%f, upper=%f, lower=%f\n", 
        //     layer, z, upper, lower);

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
  // printf("setting hint_l to %u\n", hint_l);

  const F alpha = (z - lower) / (upper - lower), 
           beta = 1.0 - alpha;

  // reset values before interpolation
  memset(v, 0, sizeof(F)*3);
  memset(f, 0, sizeof(F)*nattrs);
  if (vv)
    *vv = 0.0;

  // interpolation
  for (int i = 0; i < nverts; i ++) {
    for (int k = 0; k < 3; k ++) {
      const double Vu = V[ k + 3 * (iv[i] * nlayers + layer) ], 
	           Vl = V[ k + 3 * (iv[i] * nlayers + layer + 1) ];
      if (Vu > 1e10 || Vl > 1e10) // fill value
        return false;
      v[k] += omega[i] * ( alpha * Vu + beta * Vl );
    }

    for (int k = 0; k < nattrs; k ++)
      f[k] += omega[i] * (
                alpha * A[ k + nattrs * (iv[i] * nlayers + layer) ]
              + beta  * A[ k + nattrs * (iv[i] * nlayers + layer + 1) ]);

    if (vv) {
      *vv +=   alpha * Vv[ iv[i] * (nlayers + 1) + layer ]
             + beta  * Vv[ iv[i] * (nlayers + 1) + layer + 1];
    }
  }

  return true;
}

template <typename F=double, typename Fm=double, typename Fv=double>
__device__ __host__ 
inline static bool mpas_eval(
    const F *x,             // location
    F *v,                   // return velocity
    F *vv,                  // vertical velocity
    F *f,                   // scalar attributs
    const F alpha,          // temporal interpolation coef
    const Fv *const V[2],    // velocity field
    const Fv *const Vv[2],   // vertical velocities
    const Fv *const zTop[2], // top layer depth
    const int nattrs,             // number of scalar attributes
    const Fv *const A[2],    // scalar attributes
    const Fm *Xv,            // vertex locations
    const int max_edges,
    const int *nedges_on_cell, 
    const int *cells_on_cell,
    const int *verts_on_cell,
    const int nlayers,
    unsigned long long &hint_c, 
    unsigned int &hint_l)        // hint for searching cell and layer
{
  if (V[1]) { // two timesteps
    F _v[2][3], _vv[2], _f[2][MAX_ATTRS];

    for (int i = 0; i < 2; i ++) {
      // printf("i=%d, nattrs=%d, nlayers=%d, A[i]=%p\n", i, nattrs, nlayers, A[i]);
      if (!mpas_eval_static<F, Fm, Fv>(x, 
          _v[i], &_vv[i], _f[i], 
          V[i], Vv[i], zTop[i], nattrs, A[i], 
          Xv, max_edges, nedges_on_cell, cells_on_cell, verts_on_cell, 
          nlayers, hint_c, hint_l))
        return false;
    }

    // temporal interpolation
    const F beta = 1.0 - alpha;
    for (int j = 0; j < 3; j ++)
      v[j] = alpha * _v[0][j] + beta * _v[1][j];
    *vv = alpha * _vv[0] + beta * _vv[1];
    for (int j = 0; j < nattrs; j ++)
      f[j] = alpha * _f[0][j] + beta * _f[1][j];

    return true;
  } else
    return mpas_eval_static<F, Fm, Fv>(
        x, v, vv, f, V[0], Vv[0], zTop[0], nattrs, A[0], 
        Xv, max_edges, nedges_on_cell, cells_on_cell, verts_on_cell, 
        nlayers, hint_c, hint_l);
}

template <int ORDER=1, typename F=double, typename Fm=double, typename Fv=double>
__device__
inline static bool spherical_rk_with_vertical_velocity(
    const F h,
    const int nsteps,
    const int istep,
    ftk::feature_point_lite_t& p, 
    F *v0,                  // return velocity
    F *vv,                  // vertical velocity
    F *f,                   // scalar attributs
    const Fv *const V[2],    // velocity field
    const Fv *const Vv[2],   // vertical velocities
    const Fv *const zTop[2], // top layer depth
    const int nattrs,            // number of scalar attributes
    const Fv *const A[2],    // scalar attributes
    const Fm *const Xv,      // vertex locations
    const int max_edges,
    const int *nedges_on_cell, 
    const int *cells_on_cell,
    const int *verts_on_cell,
    const int nlayers,
    unsigned long long &hint_c, 
    unsigned int &hint_l)        // hint for searching cell and layer
{
  if (ORDER == 1) {
    const F alpha = (F)istep / nsteps;

    F v[4];
    if (!mpas_eval(p.x, v, vv, f, alpha,
          V, Vv, zTop, nattrs, A, 
          Xv, max_edges, nedges_on_cell, cells_on_cell, verts_on_cell, 
          nlayers, hint_c, hint_l))
      return false;

    for (int k = 0; k < 3; k ++)
      v0[k] = v[k];

    ftk::angular_stepping(p.x, v, h, p.x);
#if 0
    const F R = ftk::vector_2norm<3, F>(p.x); // radius
    const F R1 = R + vv[0] * h; // new radius
    for (int k = 0; k < 3; k ++)
      p.x[k] = p.x[k] * R1 / R;
#endif

    return true;
  } else if (ORDER == 4) {
    return true;
  } else
    return false;
  
  // printf("x=%f, %f, %f, v=%f, %f, %f, %f\n", 
  //     p.x[0], p.x[1], p.x[2], v[0], v[1], v[2], v[3]);
}

template <typename F=double, typename Fm=double>
__global__
static void initialize_c2v(
    const int ncells,
    const int nverts,
    const Fm *Xc, // cell x
    const Fm *Xv, // vertex x
    const int *cells_on_vert,
    F *interpolants,
    bool *vert_on_boundary)
{
  unsigned long long i = getGlobalIdx_3D_1D();
  if (i >= nverts) return;

  bool boundary = false;
  F xc[3][3], xv[3];

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
    vert_on_boundary[vi] = true;
    for (int k = 0; k < 3; k ++) // just make initchecker happy
      interpolants[k + vi*3] = 0.0; 
  } else {
    vert_on_boundary[vi] = false;
    
    F mu[3];
    bool succ = ftk::inverse_lerp_s2v3_0(xc, xv, mu);
    for (int k = 0; k < 3; k ++)
      interpolants[k + vi*3] = mu[k];
  }
}

template <typename Fv=double>
__global__
static void assign_attrs(
    Fv *A, // all attributes
    const Fv *a, // one of the attributes
    const int iattr, 
    const int nattrs,
    const int nlayers,
    const int ncells)
{
  unsigned long long i = getGlobalIdx_3D_1D();
  if (i >= ncells) return;
  const int ci = i;

  for (int layer = 0; layer < nlayers; layer ++)
    A[iattr + nattrs * (layer + ci * nlayers)] = a[layer + ci * nlayers];
}

template <typename F=double, typename Fm=double, typename Fv=double>
__global__
static void interpolate_e2c(
    Fv *Vc, // velocity on cell; please memset to zero before calling this func
    const Fv *Ve, // normal velocity on edge
    const int nlayers,
    const int nedges,
    const int ncells,
    const F *interpolants,
    const int max_edges,
    const int *nedges_on_cell,
    const int *edges_on_cell)
{
  unsigned long long i = getGlobalIdx_3D_1D();
  if (i >= ncells) return;
  
  const int ci = i;
  const int ne = nedges_on_cell[ci];

  for (int layer = 0; layer < nlayers; layer ++) {
    // for (int k = 0; k < 3; k ++) 
    //   Vc[k + 3 * (layer + ci *nlayers)] = 0.0;

    for (int j = 0; j < ne; j ++) {
      const int ei = e2ei( edges_on_cell[j + ci*max_edges] );
      // printf("%f\n", ve);
      // printf("ci=%d, ei=%d, v=%p\n", ci, ei, Ve); // layer + ei*nlayers]);
      const F ve = Ve[layer + ei * nlayers];

      for (int k = 0; k < 3; k ++) 
        Vc[k + 3 * (layer + ci *nlayers)] += ve * interpolants[k + 3 * (j + max_edges * ci)];
#if 0
      if (ci == 0)
        printf("ve=%f, vc=%f, %f, %f\n", ve, 
          Vc[3*(layer+ci*nlayers)],
          Vc[1+3*(layer+ci*nlayers)],
          Vc[2+3*(layer+ci*nlayers)]);
#endif
    }
  }
}

template <typename F=double, typename Fv=double>
__global__
static void interpolate_c2v(
    Fv *V, // vertexwise data
    const Fv *C, // cellwise data
    const int nch,
    const int nlayers,
    const int ncells,
    const int nverts,
    const F *interpolants,
    const int *cells_on_vert,
    const bool *vert_on_boundary)
{
  // unsigned long long i = 192;
  // if (threadIdx.x > 0) return;
  unsigned long long i = getGlobalIdx_3D_1D();
  if (i >= nverts) return;
  const int vi = i;

  const bool boundary = vert_on_boundary[vi];
  F lambda[3];
  int ci[3];

  for (int l = 0; l < 3; l ++) {
    lambda[l] = interpolants[l + vi * 3];
    ci[l] = c2ci(cells_on_vert[l + vi * 3]);
  }

#if 0
  if (1) { // if (blockIdx.x == 0 && blockIdx.y == 481 && threadIdx.x == 192) {
    printf("ci=%d, %d, %d, boundary=%d, lambda=%f, %f, %f\n", 
        ci[0], ci[1], ci[2], boundary, 
        lambda[0], lambda[1], lambda[2]);
  }
#endif

  for (int layer = 0; layer < nlayers; layer ++) {
    for (int k = 0; k < nch; k ++) {
      Fv &v = V[k + nch * (layer + vi * nlayers)];
      if (boundary) {
          v = 0;
      } else {
        bool invalid = false;
        v = 0;
        for (int l = 0; l < 3; l ++) {
          F val = C[k + nch * (layer + ci[l] * nlayers)];
          if (val < -1e32) {
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

template <typename F=double, typename Fm=double>
__global__
static void seed_latlonz(
      ftk::feature_point_lite_t* particles,
      const int nlat, const F lat0, const F lat1,
      const int nlon, const F lon0, const F lon1,
      const int nz,   const F z0,   const F z1, 
      const int ncells, 
      const int *kdheap,
      const Fm *Xc,
      const Fm *Xv,
      const int max_edges,
      const int *nedges_on_cell, // also n_verts_on_cell
      const int *cells_on_cell,
      const int *verts_on_cell)
{
  unsigned long long id = getGlobalIdx_3D_1D();
  if (id >= nlat * nlon * nz) return;
  
  const F dz   = nz == 1 ? 0.0 : (z1 - z0) / (nz - 1), 
          dlat = nlat == 1 ? 0.0 : (lat1 - lat0) / (nlat - 1),
          dlon = nlon == 1 ? 0.0 : (lon1 - lon0) / (nlon - 1);

  const int k = id % nz;
  const int j = (id / nz) % nlon;
  const int i = id / (nlon*nz);
        
  const F lat = ftk::deg2rad(i * dlat + lat0);
  const F slat = sin(lat), 
          clat = cos(lat);

  const F lon = ftk::deg2rad(j * dlon + lon0);
  const F clon = cos(lon),
          slon = sin(lon);

  const F z = k * dz + z0;
  const F r = R0 + z;

  ftk::feature_point_lite_t &p = particles[id];
  p.x[0] = r * clon * clat;
  p.x[1] = r * slon * clat;
  p.x[2] = r * slat;
  p.t = 0.0;
  
  int iv[7];
  F xv[7][3];
  p.tag /* hint_c */ = locate_cell_global(p.x, iv, xv, ncells, kdheap, 
      Xc, Xv, max_edges, nedges_on_cell, cells_on_cell, verts_on_cell);
}

template <typename F=double, typename Fm=double, typename Fv=double>
__global__
static void mpas_trace(
    const F h,
    const int nsteps,
    const int nsubsteps,
    const int isubstart,
    const int nparticles,
    ftk::feature_point_lite_t* particles,
    const Fv *const V[2],    // velocity field
    const Fv *const Vv[2],   // vertical velocities
    const Fv *const zTop[2], // top layer depth
    const int nattrs,   // number of scalar attributes
    const Fv *const A[2],    // scalar attributes
    const Fm *Xv,   // vertex locations
    const int max_edges,
    const int *nedges_on_cell, 
    const int *cells_on_cell,
    const int *verts_on_cell,
    const int nlayers)
{
  unsigned long long i = getGlobalIdx_3D_1D();
  if (i >= nparticles) return;

#if 0 // debuggin a specifc particle
  if (getGlobalIdx_3D_1D() > 0) return;
  int bid = 155;
  int tid = bid * 167 + 130;
  const auto i = tid; 
  //(155,0,0), thread (130,0,0) in (167,1,1),(256,1,1);
#endif

  ftk::feature_point_lite_t &p = particles[i];
  F v0[3], vv[3], f[MAX_ATTRS];
  // const D h = 1.0 / nsteps;

  unsigned long long &hint_c = p.tag;
  unsigned int &hint_l = p.type;

  for (int j = isubstart; j < nsubsteps; j ++) {
    bool succ = spherical_rk_with_vertical_velocity<1>(
        h, nsteps, j, 
        p, v0, vv, f, 
        V, Vv, zTop, nattrs, A, Xv, 
        max_edges, nedges_on_cell, cells_on_cell, verts_on_cell, 
        nlayers, hint_c, hint_l);

    if (succ) {
      for (int k = 0; k < 3; k ++) {
        p.scalar[k] = v0[k];
      }
      p.scalar[3] = vv[0];
      for (int k = 0; k < nattrs; k ++) {
        p.scalar[k+4] = f[k];
      }
    } else
      break;
  }
}

///////////////////////////
void mop_create_ctx(mop_ctx_t **c_, int device,
    bool prec_compute, bool prec_mesh, bool prec_var)
{
  *c_ = (ctx_t*)malloc(sizeof(ctx_t));
  ctx_t *c = *c_;
  memset(c, 0, sizeof(ctx_t));

  c->prec_compute = prec_compute;
  c->prec_mesh = prec_mesh;
  c->prec_var = prec_var;

  fprintf(stderr, "prec: compute (%d), mesh (%d), vars (%d)\n", 
      prec_compute, prec_mesh, prec_var);

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
  if (c->d_cells_on_edge != NULL) cudaFree(c->d_cells_on_edge);
  if (c->d_cells_on_vert != NULL) cudaFree(c->d_cells_on_vert);
  if (c->d_edges_on_cell != NULL) cudaFree(c->d_edges_on_cell);
  if (c->d_verts_on_cell != NULL) cudaFree(c->d_verts_on_cell);

  if (c->d_kdheap != NULL) cudaFree(c->d_kdheap);

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
  if (c->d_vert_on_boundary != NULL) cudaFree(c->d_vert_on_boundary);
  if (c->dcw != NULL) cudaFree(c->dcw);
  if (c->dew != NULL) cudaFree(c->dew);

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
    *hbuf = (T*)realloc(*hbuf, m * sizeof(T));

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
 
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] load particles");
}

void mop_load_mesh(mop_ctx_t *c,
    const int ncells,
    const int nedges,
    const int nverts, 
    // const int nlayers, 
    const int max_edges,
    const int nattrs,
    const void *Xc,
    const void *Xv,
    const int *nedges_on_cell, 
    const int *cells_on_cell,
    const int *cells_on_edge,
    const int *cells_on_vert,
    const int *edges_on_cell,
    const int *verts_on_cell)
{
  c->ncells = ncells;
  // c->nlayers = nlayers;
  c->nedges = nedges;
  c->nverts = nverts;
  c->max_edges = max_edges;
  c->nattrs = nattrs;

  size_t fsize = c->prec_mesh ? sizeof(double) : sizeof(float);

  cudaMalloc((void**)&c->d_Xc, 
      size_t(ncells) * fsize * 3);
  cudaMemcpy(c->d_Xc, Xc, 
      size_t(ncells) * fsize * 3,
      cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_Xv, 
      size_t(nverts) * fsize * 3);
  cudaMemcpy(c->d_Xv, Xv, 
      size_t(nverts) * fsize * 3, 
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
  
  cudaMalloc((void**)&c->d_cells_on_edge, 
      size_t(nedges) * 2 * sizeof(int));
  cudaMemcpy(c->d_cells_on_edge, cells_on_edge, 
      size_t(nedges) * 2 * sizeof(int),
      cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_cells_on_vert, 
      size_t(nverts) * 3 * sizeof(int));
  cudaMemcpy(c->d_cells_on_vert, cells_on_vert, 
      size_t(nverts) * 3 * sizeof(int),
      cudaMemcpyHostToDevice);
  
  cudaMalloc((void**)&c->d_edges_on_cell, 
      size_t(ncells) * max_edges * sizeof(int));
  cudaMemcpy(c->d_edges_on_cell, edges_on_cell, 
      size_t(ncells) * max_edges * sizeof(int), 
      cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_verts_on_cell, 
      size_t(ncells) * max_edges * sizeof(int));
  cudaMemcpy(c->d_verts_on_cell, verts_on_cell, 
      size_t(ncells) * max_edges * sizeof(int), 
      cudaMemcpyHostToDevice);

  checkLastCudaError("[FTK-CUDA] loading mpas mesh");

  // initialize interpolants
  {
    cudaMalloc((void**)&c->d_c2v_interpolants, 
        size_t(c->nverts) * 3 * sizeof(double)); // interpolants are in double, at least for now

    cudaMalloc((void**)&c->d_vert_on_boundary,
        size_t(c->nverts) * sizeof(bool));
    
    // cudaDeviceSynchronize();
    checkLastCudaError("[FTK-CUDA] initializing mpas c2v interpolants: malloc");

    size_t ntasks = c->nverts;
    const int nBlocks = idivup(ntasks, blockSize);
    dim3 gridSize;
  
    if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
    else gridSize = dim3(nBlocks);

    fprintf(stderr, "initializing c2v: ncells=%d, nverts=%d\n",
        c->ncells, c->nverts);

    if (c->prec_mesh) {
      initialize_c2v<double, double><<<gridSize, blockSize>>>(
          c->ncells, 
          c->nverts, 
          (double*)c->d_Xc,
          (double*)c->d_Xv,
          c->d_cells_on_vert, 
          (double*)c->d_c2v_interpolants,
          c->d_vert_on_boundary);
    } else {
      initialize_c2v<double, float><<<gridSize, blockSize>>>(
          c->ncells, 
          c->nverts, 
          (float*)c->d_Xc,
          (float*)c->d_Xv,
          c->d_cells_on_vert, 
          (double*)c->d_c2v_interpolants,
          c->d_vert_on_boundary);
    }
 
    fprintf(stderr, "c2v initialization done.\n");
    // cudaDeviceSynchronize();
    checkLastCudaError("[FTK-CUDA] initializing mpas c2v interpolants");
  }

  // initialize normal vectors
  {

  }
}

void mop_load_e2c_interpolants(mop_ctx_t *c, const void *p)
{
  size_t fsize = sizeof(double); 

  cudaMalloc((void**)&c->d_e2c_interpolants,
      size_t(c->ncells) * c->max_edges * fsize * 3);

  cudaMemcpy(c->d_e2c_interpolants, p,
      size_t(c->ncells) * c->max_edges * fsize * 3, cudaMemcpyHostToDevice);

  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] load e2c interpolants");
}

static void load_data(
    void** dbuf,
    const void *buf,
    const size_t fsize,
    const size_t n, 
    const char *name)
{
  void *d;
  if (dbuf[0] == NULL) {
    cudaMalloc((void**)&dbuf[0], fsize * n);
    checkLastCudaError("[FTK-CUDA] loading data 0");
    d = dbuf[0];
  } else if (dbuf[1] == NULL) {
    cudaMalloc((void**)&dbuf[1], fsize * n);
    checkLastCudaError("[FTK-CUDA] loading data 1");
    d = dbuf[1];
  } else {
    std::swap(dbuf[0], dbuf[1]);
    d = dbuf[1];
  }
  cudaMemcpy(d, buf, fsize * n, cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] cudaMemcpy");
}

static void load_cw_data(
    mop_ctx_t *c,
    void **dbuf, // an array with two device pointers
    void ***ddbuf, // device array of pointers
    const void *cw,
    bool cw_on_gpu,
    const int nch,
    const int nlayers) // host pointer with cw data
{
  size_t fsize = c->prec_var ? sizeof(double) : sizeof(float);
  
  // initialize buffers for vertexwise data in two adjacent timesteps
  // cudaDeviceSynchronize();
  checkLastCudaError("memcpy to c2w buffer: load cw 0");

  void *d;
  if (dbuf[0] == NULL) {
    cudaMalloc((void**)&dbuf[0], 
        fsize * c->nverts * nch * nlayers);
    d = dbuf[0];
    checkLastCudaError("malloc vw buffers 0");
  } else if (dbuf[1] == NULL) {
    cudaMalloc((void**)&dbuf[1], 
        fsize * c->nverts * nch * nlayers);
    d = dbuf[1];
    checkLastCudaError("malloc vw buffers 1");
  } else {
    std::swap(dbuf[0], dbuf[1]);
    d = dbuf[1];
  }
 
  // cudaDeviceSynchronize();
  checkLastCudaError("memcpy to c2w buffer: dev ptrs0");
  if (*ddbuf == NULL) {
    cudaMalloc((void**)ddbuf, sizeof(void*) * 2);
    checkLastCudaError("memcpy to c2w buffer: allocating ddbuf");
  }
  // fprintf(stderr, "ddbuf=%p, dbuf=%p\n", *ddbuf, dbuf);
  cudaMemcpy(*ddbuf, dbuf, sizeof(void*) * 2, 
      cudaMemcpyHostToDevice);
  checkLastCudaError("memcpy to c2w buffer: dev ptrs");
  
  // copy data to c2w buffer
  void const* dcw; 
  if (cw_on_gpu) {
    dcw = cw; // cw data is already on gpu
  } else {
    dcw = c->dcw;
    assert(c->dcw != NULL);
    cudaMemcpy(c->dcw, cw, 
        fsize * c->ncells * nch * nlayers, 
        cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    // fprintf(stderr, "dcw=%p\n", c->dcw);
    checkLastCudaError("memcpy to c2w buffer");
  }

  // c2w interpolation
  {
    size_t ntasks = c->nverts;
    const int nBlocks = idivup(ntasks, blockSize);
    dim3 gridSize;
  
    if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
    else gridSize = dim3(nBlocks);

    // cudaDeviceSynchronize();
    if (c->prec_var) { 
      interpolate_c2v<double, double><<<gridSize, blockSize>>>(
          (double*)d, (double*)dcw, nch, nlayers,
          c->ncells,
          c->nverts,
          (double*)c->d_c2v_interpolants,
          c->d_cells_on_vert,
          c->d_vert_on_boundary);
    } else {
      interpolate_c2v<double, float><<<gridSize, blockSize>>>(
          (float*)d, (float*)dcw, nch, nlayers,
          c->ncells,
          c->nverts,
          (double*)c->d_c2v_interpolants,
          c->d_cells_on_vert,
          c->d_vert_on_boundary);
    }
 
    // cudaDeviceSynchronize();
    checkLastCudaError("c2w interpolation");
  }
}

void mop_load_kdheap(mop_ctx_t *c, int *heap)
{
  if (c->d_kdheap == NULL)
    cudaMalloc((void**)&c->d_kdheap, c->ncells * sizeof(int));

  cudaMemcpy(c->d_kdheap, heap, c->ncells * sizeof(int), 
      cudaMemcpyHostToDevice);
  
  checkLastCudaError("load kdheap");
}

void mop_load_data(mop_ctx_t *c, 
    const double *V,
    const double *Vv,
    const double *zTop,
    const double *A)
{
  size_t fsize = c->prec_var ? sizeof(double) : sizeof(float);
  
  load_data(c->d_V, V, fsize, 3 * c->nverts * c->nlayers, "V");
  load_data(c->d_Vv, Vv, fsize, c->nverts * (c->nlayers+1), "Vv");
  load_data(c->d_zTop, zTop, fsize, c->nverts * c->nlayers, "zTop");
  load_data(c->d_A, A, fsize, c->nattrs * c->nverts * c->nlayers, "f");
}

void mop_initialize_dcw(mop_ctx_t *c)
{
  size_t fsize = c->prec_var ? sizeof(double) : sizeof(float);
  
  if (c->dcw == NULL)
    cudaMalloc((void**)&c->dcw, 
        fsize * c->ncells * std::max(3, c->nattrs) * (c->nlayers+1)); // a sufficiently large buffer for loading cellwise data
  // cudaDeviceSynchronize();
  checkLastCudaError("malloc dcw");
}

void mop_initialize_dew(mop_ctx_t *c)
{
  size_t fsize = c->prec_var ? sizeof(double) : sizeof(float);
  
  if (c->dew == NULL)
    cudaMalloc((void**)&c->dew, 
        fsize * std::max(c->nedges, c->ncells) * std::max(3, c->nattrs) * c->nlayers); // a sufficiently large buffer for loading edgewise data
        // sizeof(double) * c->nedges * c->nlayers); // a sufficiently large buffer for loading edgewise data
  // cudaDeviceSynchronize();
  checkLastCudaError("malloc dew");
}

void mop_load_data_cw(mop_ctx_t *c,
    const double t,
    const void *Vc,
    const void *Vvc,
    const void *zTopc,
    const void *Ac)
{
  mop_initialize_dcw(c);

  std::swap(c->T[0], c->T[1]);
  c->T[0] = t;

  load_cw_data(c, c->d_V, &c->dd_V, Vc, false, 3, c->nlayers); // V
  load_cw_data(c, c->d_Vv, &c->dd_Vv, Vvc, false, 1, c->nlayers+1); // Vv
  load_cw_data(c, c->d_zTop, &c->dd_zTop, zTopc, false, 1, c->nlayers); // zTop
  load_cw_data(c, c->d_A, &c->dd_A, Ac, false, 2 /*c->nattrs*/, c->nlayers); // f
}

void mop_load_data_with_normal_velocity(mop_ctx_t *c,
    const double t,
    const void *V, // normal velocity
    const void *Vvc, 
    const void *zTopc,
    const void *Ac[])
{
  size_t fsize = c->prec_var ? sizeof(double) : sizeof(float);
  
  mop_initialize_dcw(c);
  mop_initialize_dew(c);

  std::swap(c->T[0], c->T[1]);
  c->T[0] = t;
    
  size_t ntasks = c->ncells;
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;

  if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else gridSize = dim3(nBlocks);

  // interpolate e2c
  {
    // cudaDeviceSynchronize();
    checkLastCudaError("e2c interpolation: 0");
    cudaMemcpy(c->dew, V, fsize * c->nedges * c->nlayers, cudaMemcpyHostToDevice);
    cudaMemset(c->dcw, 0, fsize * 3 * c->ncells * c->nlayers);
    
    // cudaDeviceSynchronize();
    checkLastCudaError("e2c interpolation: 1");
   
    if (c->prec_var) {
      interpolate_e2c<double, double, double><<<gridSize, blockSize>>>(
         (double*)c->dcw, (double*)c->dew, c->nlayers, c->nedges, c->ncells,
         (double*)c->d_e2c_interpolants, c->max_edges, 
         c->d_nedges_on_cell, c->d_edges_on_cell);
    } else {
      interpolate_e2c<double, double, float><<<gridSize, blockSize>>>(
         (float*)c->dcw, (float*)c->dew, c->nlayers, c->nedges, c->ncells,
         (double*)c->d_e2c_interpolants, c->max_edges, 
         c->d_nedges_on_cell, c->d_edges_on_cell);
    }
    
    // cudaDeviceSynchronize();
    checkLastCudaError("e2c interpolation");
  }

  // fprintf(stderr, "loading V..\n");
  load_cw_data(c, c->d_V, &c->dd_V, c->dcw, true, 3, c->nlayers); // V
  // fprintf(stderr, "loading Vv..\n");
  load_cw_data(c, c->d_Vv, &c->dd_Vv, Vvc, false, 1, c->nlayers+1); // Vv
  // fprintf(stderr, "loading zTop..\n");
  load_cw_data(c, c->d_zTop, &c->dd_zTop, zTopc, false, 1, c->nlayers); // zTop
 
  void *tmp = c->dew;
  fprintf(stderr, "loading attrs..\n");
  for (int i = 0; i < c->nattrs; i ++) {
    cudaMemcpy(tmp, Ac[i], fsize * c->ncells * c->nlayers, cudaMemcpyHostToDevice);

    if (c->prec_var) 
      assign_attrs<double><<<gridSize, blockSize>>>((double*)c->dcw, (double*)tmp, i, c->nattrs, c->nlayers, c->ncells);
    else
      assign_attrs<float><<<gridSize, blockSize>>>((float*)c->dcw, (float*)tmp, i, c->nattrs, c->nlayers, c->ncells);
    
    checkLastCudaError("load attrs: assign");
  }
  load_cw_data(c, c->d_A, &c->dd_A, c->dcw, true, c->nattrs, c->nlayers); // f
}

void mop_seed_latlonz(mop_ctx_t *c,
      const int nlat, const double lat0, const double lat1,
      const int nlon, const double lon0, const double lon1,
      const int nz,   const double z0,   const double z1)
{
  size_t ntasks = nlat * nlon * nz;
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;

  if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else gridSize = dim3(nBlocks);
  
  // allocate
  realloc_both(&c->hparts, &c->dparts, c->nparticles, ntasks);
  c->nparticles = ntasks;

  // seed
  seed_latlonz<double, double><<<gridSize, blockSize>>>(
      c->dparts, 
      nlat, lat0, lat1, 
      nlon, lon0, lon1,
      nz, z0, z1,
      c->ncells,
      c->d_kdheap,
      (double*)c->d_Xc, 
      (double*)c->d_Xv, 
      c->max_edges,
      c->d_nedges_on_cell,
      c->d_cells_on_cell,
      c->d_verts_on_cell);
  
  cudaMemcpy(c->hparts, c->dparts, 
      c->nparticles * sizeof(ftk::feature_point_lite_t), 
      cudaMemcpyDeviceToHost);
  
  checkLastCudaError("[FTK-CUDA] mop_seed_latlonz");
}

void mop_execute(mop_ctx_t *c,
    const double T, // duration
    const int nsteps, 
    const int nsubsteps)
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
  // fprintf(stderr, "%p, %p\n", c->d_V[0], c->d_V[1]);

#if 0
  const double duration = 86400.0 * 365 * 10; // 10 years
  const int nsteps = 1024 * 365 * 10; // 1024; // 32768;
  const int nsubsteps = nsteps;
  const double h = duration / nsteps;
#endif
  
  const double h = T / nsteps;

  // fprintf(stderr, "h=%f, %p, %p, %p, %p, %d, %p, %p, %p, %p, %p\n", 
  //     h, c->dparts, c->dd_V, c->dd_Vv, c->dd_zTop, c->nattrs,
  //     c->dd_A, c->d_Xv, c->d_nedges_on_cell, c->d_cells_on_cell, c->d_verts_on_cell);

  if (c->prec_mesh) {
    if (c->prec_var) {
      mpas_trace<double, double, double><<<gridSize, blockSize>>>(
          h,
          nsteps, nsubsteps, 0,
          ntasks, 
          c->dparts,
          (double**)c->dd_V,
          (double**)c->dd_Vv,
          (double**)c->dd_zTop, 
          c->nattrs, 
          (double**)c->dd_A,
          (double*)c->d_Xv,
          c->max_edges, 
          c->d_nedges_on_cell, 
          c->d_cells_on_cell,
          c->d_verts_on_cell, 
          c->nlayers);
    } else {
      mpas_trace<double, double, float><<<gridSize, blockSize>>>(
          h,
          nsteps, nsubsteps, 0,
          ntasks, 
          c->dparts,
          (float**)c->dd_V,
          (float**)c->dd_Vv,
          (float**)c->dd_zTop, 
          c->nattrs, 
          (float**)c->dd_A,
          (double*)c->d_Xv,
          c->max_edges, 
          c->d_nedges_on_cell, 
          c->d_cells_on_cell,
          c->d_verts_on_cell, 
          c->nlayers);
    }
  } else {
    if (c->prec_var) assert(false);
    else {
      mpas_trace<double, float, float><<<gridSize, blockSize>>>(
          h,
          nsteps, nsubsteps, 0,
          ntasks, 
          c->dparts,
          (float**)c->dd_V,
          (float**)c->dd_Vv,
          (float**)c->dd_zTop, 
          c->nattrs, 
          (float**)c->dd_A,
          (float*)c->d_Xv,
          c->max_edges, 
          c->d_nedges_on_cell, 
          c->d_cells_on_cell,
          c->d_verts_on_cell, 
          c->nlayers);
    }
  }
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] mop_execute");

  cudaMemcpy(c->hparts, c->dparts, 
      c->nparticles * sizeof(ftk::feature_point_lite_t), 
      cudaMemcpyDeviceToHost);
  checkLastCudaError("[FTK-CUDA] mop_execute: memcpy");

  // fprintf(stderr, "exiting kernel\n");
}

void mop_swap(mop_ctx_t *c)
{
  // std::swap(c->d_V[0], c->d_V[1]); // make no sense for now
}
