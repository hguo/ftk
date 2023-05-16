#include <nvfunctional>
#include <ftk/numeric/mpas.hh>
#include <ftk/numeric/wachspress_interpolation.hh>

static const int MAX_VERTS = 10;
static const int MAX_LAYERS = 100;
static const double R0 = 6371229.0;

// what we need:
// - velocity field
// - vertical velocity field
// - zTop field
// - scalar attribute fields

template <typename F=double>
__device__ __host__
inline bool point_in_mpas_cell(
    const int cell,
    const double *x,
    int iv[],
    double xv[][3], // returns vertex coordinates
    const int max_edges,
    const double *Xv,
    const int *n_edges_on_cell, 
    const int *verts_on_cell)
{
  // if (cell < 0) return false;
  const nverts = n_edges_on_cell[cell];
  // double xv[MAX_VERTS][3];

  for (int i = 0; i < nverts; i ++) {
    for (int k = 0; k < 3; k ++) {
      const int vertex = verts_on_cell[cell * nverts + i];
      iv[i][k] = vertex;
      xv[i][k] = Xv[vertex*3+k];
    }
  }

  return point_in_mpas_cell(nverts, xv, x);
}

__device__ __host__
static int locate_cell_local( // local search among neighbors
    const int curr, // current cell
    const double *x,
    double iv[], // returns vertex ids
    double xv[][3], // returns vertex coordinates
    const double *Xv,
    const int max_edges,
    const int *n_edges_on_cell, // also n_verts_on_cell
    const int *cells_on_cell,
    const int *verts_on_cell)
{
  if (curr < 0)
    return -1; // not found
  else if (point_in_mpas_cell(
        curr, x, iv, xv, 
        max_edges, Xv, n_edges_on_cell, verts_on_cell))
    return curr;
  else {
    for (int i = 0; i < n_edges_on_cell[curr]; i ++) {
      const int cell = cells_on_cell[i + max_edges * curr];
      if (point_in_mpas_cell(
            cell, x, iv, xv,
            max_edges, Xv, n_edges_on_cell, verts_on_cell))
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
    const int nattrs,   // number of scalar attributes
    const double *A,    // scalar attributes
    const double *Xv,   // vertex locations
    const int max_edges,
    const int *n_edges_on_cell, 
    const int *cells_on_cell,
    const int *verts_on_cell,
    const int nlayers,
    const double *zTop, // top layer depth
    int &hint_c, 
    int &hint_l)        // hint for searching cell and layer
{
  int iv[MAX_VERTS];
  double xv[MAX_VERTS][3];

  const int cell = locate_cell_local(hint_c, 
      x, iv, xv, 
      Xv, max_edges, n_edges_on_cell, 
      cells_on_cell, verts_on_cell);
  if (cell < 0) return false;
  else hint_c = cell;

  const int nverts = n_edges_on_cell[cell];

  // compute weights based on xyzVerts
  double omega[MAX_VERTS]; 
  wachspress_weights(nverts, xv, x, omega); 

  // locate layer
  int layer = hint_l;
  double upper = 0.0, lower = 0.0;
  const double R = vector_2norm<3>(x);
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

  if (!succ) {
    if (dir == 1) { // downward
      upper = lower;
      for (layer = layer + 1 ; layer < nlayers-1; layer ++) {
        lower = 0.0;
        for (int k = 0; k < nverts; k ++)
          lower += omega[k] * zTop[ iv[i] * nlayers + layer + 1];

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
          upper += omega[k] * zTop[ iv[i] * nlayers + layer];

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
               beta = 1.0 = alpha;

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

    for (int k = 0; k < nattrs; k ++)
      f[k] += omega[i] * (
                alpha * A[ k + nattrs * (iv[i] * nlayers + layer) ]
              + beta  * A[ k + nattrs * (iv[i] * nlayers + layer + 1) ]);

    *vv +=   alpha * Vv[ iv[i] * (nlayers + 1) + layer ]
           + beta  * Vv[ iv[i] * (nlayers + 1) + layer + 1];
  }

  return true;
}
