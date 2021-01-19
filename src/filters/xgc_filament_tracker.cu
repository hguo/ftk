// #include <nvfunctional>
// #include "common.cuh"
#include "mx4.cuh"

// what are needed in device memory:
// - field-following interpolants for virtual poloidal planes
// - triangles, edges, and vertex coordinates in 2D mesh
// - current and next timestep of scalar, vector, and jacobian fields

template <typename I, typename F>
__device__
bool check_simplex(
    I current_timestep, 
    I i,
    const I nphi,
    const I iphi,
    const I vphi,
    const I m2n0,
    const I m2n1,
    const I m2n2,
    const F m2coords[], // vertex coordinates
    const I m2edges[], 
    const I m2tris[],
    const xgc_interpolant_t **interpolants, // interpolants
    const F *scalar[2], // current and next scalar
    const F *vector[2], // current and next grad
    const F *jacobian[2]) // current and next jacobian
    // cp_t & cp) // WIP: critical points
{
  // typedef ftk::fixed_point<> fp_t;
  typedef unsigned long long fp_t; // WIP

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

  fp_t vf[3][2];
}
