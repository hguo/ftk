#ifndef _FTK_MX3_CUH
#define _FTK_MX3_CUH

#include "me2.cuh"

struct xgc_interpolant_t {
  int tri0[3], tri1[3];
  double mu0[3], mu1[3];
};

template <typename I>
__device__
I mx3_transform_vert(I i, const I mx3n0)
{
  while (i < 0) i += mx3n0;
  return i % mx3n0;
}

template <typename I, typename F>
__device__
inline void mx3_get_coords(
    I v3, F x[], 
    const I m2n0, const F m2coords[])
{
  return me2_get_coords(v3, x, m2n0, m2coords);
}

template <typename I>
__device__
void mx3_get_edge(I k, I verts[2], 
    const I np, const I m2n0, const I m2n1, 
    const I m2edges[])
{
  const I mx3n0 = m2n0 * np;
  me2_get_edge(k, verts, m2n0, m2n1, m2edges);
  for (int i = 0; i < 2; i ++)
    verts[i] = mx3_transform_vert( verts[i], mx3n0 );

  sort2(verts);
}

template <typename I>
__device__
void mx3_get_tri(I k, I verts[3], const 
    I np, 
    const I m2n0, const I m2n1, const I m2n2, 
    const I m2edges[], const I m2tris[])
{
  const I mx3n0 = m2n0 * np;
  me2_get_tri(k, verts, m2n0, m2n1, m2n2, m2edges, m2tris);
  for (int i = 0; i < 3; i ++) 
    verts[i] = mx3_transform_vert( verts[i], mx3n0 );

  sort3(verts);
}

template <typename I, typename F>
__device__
inline F mx3_interpolate_scalar(const I i,
    const I nphi, const I iphi, const I vphi,
    const I m2n0, 
    const ftk::xgc_interpolant_t<I, F> *interpolants,
    const F *scalar)
{
  const I p = i / m2n0; // poloidal plane
  if (p % vphi == 0) {
    const I p0 = p / vphi;
    const I k = i % m2n0;
    const I idx = m2n0 * p0 + k;

    return scalar[idx];
  } else { // virtual plane
    const auto& l = interpolants[(p % vphi - 1) * m2n0 + i % m2n0];
    if (l.tri0[0] < 0 || l.tri1[0] < 0) 
      return nan("");

    const I p0 = p / vphi, p1 = (p0 + 1) % nphi;
    const F beta = F(p) / vphi - p0, alpha = F(1) - beta;

    F f = 0;
    for (int k = 0; k < 3; k ++) {
      const I idx0 = m2n0 * p0 + l.tri0[k], 
              idx1 = m2n0 * p1 + l.tri1[k];

      f += alpha * l.mu0[k] * scalar[idx0] + beta * l.mu1[k] * scalar[idx1];
    }

    return f;
  }
}

template <typename I, typename F>
__device__
inline void mx3_interpolate(const I i,
    const I nphi, const I iphi, const I vphi,
    const I m2n0, 
    const ftk::xgc_interpolant_t<I, F> *interpolants, 
    const F *scalar, const F *vector, const F *jacobian, 
    F &f, F v[2], F j[2][2])
{
  const I p = i / m2n0; // poloidal plane
  if (p % vphi == 0) {
    const I p0 = p / vphi;
    const I k = i % m2n0;
    const I idx = m2n0 * p0 + k;

    f = scalar[idx];
    v[0] = vector[2*idx];
    v[1] = vector[2*idx+1];
    j[0][0] = jacobian[4*idx];
    j[0][1] = jacobian[4*idx+1];
    j[1][0] = jacobian[4*idx+2];
    j[1][1] = jacobian[4*idx+3];
  } else { // virtual plane
    const I p0 = p / vphi, p1 = (p0 + 1) % nphi;
    const F beta = F(p) / vphi - p0, alpha = F(1) - beta;
    // const auto& l = interpolants[p % vphi][i % m2n0];
    const auto& l = interpolants[(p % vphi - 1) * m2n0 + i % m2n0];

    // init
    f = 0;
    for (int k0 = 0; k0 < 2; k0 ++) {
      v[k0] = 0;
      for (int k1 = 0; k1 < 2; k1 ++) 
        j[k0][k1] = 0;
    }

    // add values
    for (int k = 0; k < 3; k ++) {
      const I idx0 = m2n0 * p0 + l.tri0[k], 
              idx1 = m2n0 * p1 + l.tri1[k];

      f += alpha * l.mu0[k] * scalar[idx0] + beta * l.mu1[k] * scalar[idx1];
      for (int k0 = 0; k0 < 2; k0 ++) {
        v[k0] += alpha * l.mu0[k] * vector[idx0*2+k0] + beta * l.mu1[k] * vector[idx1*2+k0];
        for (int k1 = 0; k1 < 2; k1 ++) 
          j[k0][k1] += alpha * l.mu0[k] * jacobian[idx0*4+k0*2+k1] + beta * l.mu1[k] * jacobian[idx1*4+k0*2+k1];
      }
    }
  }
}

#endif
