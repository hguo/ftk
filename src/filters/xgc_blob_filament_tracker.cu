#include <nvfunctional>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/clamp.hh>
#include <ftk/numeric/symmetric_matrix.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/filters/xgc_blob_filament_tracker.cuh>
#include "common.cuh"
#include "locator2.cuh"
#include "mx4.cuh"

// what are needed in device memory:
// - field-following interpolants for virtual poloidal planes
// - triangles, edges, and vertex coordinates in 2D mesh
// - current and next timestep of scalar, vector, and jacobian fields

using namespace ftk;

typedef xft_ctx_t ctx_t;

template <typename I, typename F>
__global__ void m2_cellwise_scalar_gradient(
    const I m2n2, 
    const I m2tris[], 
    const F m2coords[],
    const F scalar[], // input
    F cwgrad[]) // output
{
  int tid = getGlobalIdx_3D_1D();
  I i = tid;
  if (i >= m2n2) return;

  I tri[3];
  m2_get_tri(i, tri, m2tris);

  F X[3][2], f[3], gradf[2];
  for (int j = 0; j < 3; j ++) {
    const I k = tri[j];
    X[j][0] = m2coords[2*k];
    X[j][1] = m2coords[2*k+1];
    f[j] = scalar[k];
  }
  
  gradient_2dsimplex2(X, f, gradf);
  cwgrad[2*i] = gradf[0];
  cwgrad[2*i+1] = gradf[1];
}

template <typename I, typename F>
__global__ void m2_compute_invdet(
    const I m2n2, 
    const I m2tris[], 
    const F m2coords[],
    F invdet[]) // output
{
  int tid = getGlobalIdx_3D_1D();
  I i = tid;
  if (i >= m2n2) return;
  
  I tri[3];
  m2_get_tri(i, tri, m2tris);
  
  F X[3][2];
  for (int j = 0; j < 3; j ++) {
    const I k = tri[j];
    X[j][0] = m2coords[2*k];
    X[j][1] = m2coords[2*k+1];
  }

  const F det = (X[1][1] - X[2][1]) * (X[0][0] - X[2][0]) 
    + (X[2][0] - X[1][0]) * (X[0][1] - X[2][1]);

  invdet[i] = F(1) / det;
}

template <typename I, typename F>
__global__ void m2_cellwise2vertexwise_scalar_gradient(
    const I m2n0, const I max_vertex_triangles,
    const I vertex_triangles[],
    const F cwgrad[], // input
    F grad[]) // output
{
  int tid = getGlobalIdx_3D_1D();
  I i = tid;
  if (i >= m2n0) return;

  F gradf[2] = {0};
  int ntris = 0;
  for (int j = 0; j < max_vertex_triangles; j ++) {
    const I tri = vertex_triangles[i*max_vertex_triangles+j];
    if (tri >= 0) {
      gradf[0] += cwgrad[tri*2];
      gradf[1] += cwgrad[tri*2+1];
      ntris ++;
    }
  }
  grad[i*2] = gradf[0] / ntris;
  grad[i*2+1] = gradf[1] / ntris;
}

template <typename I, typename F>
__global__ void mx3_upsample_scalar(
    const I nphi, const I iphi, const I vphi,
    const I m2n0, 
    const ftk::xgc_interpolant_t<I, F>* interpolants,
    const F *input, 
    F *output)
{
  int tid = getGlobalIdx_3D_1D();
  I i = tid;
  if (i >= m2n0 * nphi * vphi) return;

  output[i] = mx3_interpolate_scalar(i, 
      nphi, iphi, vphi, m2n0, 
      interpolants, input);
}

template <typename I, typename F>
__global__ void mx3_derive_deltaB(
    const I nphi, const I iphi, const I vphi, 
    const I m2n0, 
    const F *apars, // upsampled
    const F *gradAs, // single plane
    const F *bfield,
    const F *bfield0, 
    const F *curl_bfield0,
    const I p,
    F *deltaB) // output for one single poloidal plane
{
  int tid = getGlobalIdx_3D_1D();
  I i = tid;
  if (i >= m2n0) return; 

  const I np = nphi * vphi;
  const I nip = np * iphi; // considering iphi
  const I ii = p * m2n0 + i;
  // const I p = ii / m2n0;
  // const I i = ii % m2n0;
  const I pnext = (p + 1) % np, 
          pprev = (p + np - 1) % np;
  const F dphi = 2 * M_PI / nip;

  const F dAsdphi = (apars[pnext*m2n0+i] - apars[pprev*m2n0+i]) / (dphi * 2); // central difference

  deltaB[i*3] = gradAs[i*2+1] * bfield0[i*3+2] - dAsdphi * bfield0[i*3+1] + apars[ii] * curl_bfield0[i*3];
  deltaB[i*3+1] = -gradAs[i*2] * bfield0[i*3+2] + dAsdphi * bfield0[i*3] + apars[ii] * curl_bfield0[i*3+1];
  deltaB[i*3+2] = gradAs[i*2] * bfield0[i*3+1] - gradAs[i*2+1] * bfield0[i*3] + apars[ii] * curl_bfield0[i*3+2];
}

template <typename I, typename F>
__device__
inline bool eval_static_b(
    const I m2n0,
    const I nphi, // const I iphi, const I vphi,
    const I *m2tris,
    const F *m2invdet, 
    const F *staticB,
    const bvh2d_node_t<I, F> *bvh, 
    const F rz[2], // input rz
    F v[2]) // output normalized vector
{
  F mu[3]; // barycentric coordinates
  const I tid = bvh2_locate_point_tri(bvh, rz[0], rz[1], mu, m2invdet);

  if (tid >= 0) {
    F B[3][3], b[3];
    I tri[3];
    m2_get_tri(tid, tri, m2tris);
    
    for (int i = 0; i < 3; i ++)
      for (int j = 0; j < 3; j ++)
        B[i][j] = staticB[tri[i]*3 + j];

    lerp_s2v3(B, mu, b);
    v[0] = rz[0] * b[0] / b[2];
    v[1] = rz[0] * b[1] / b[2];
    
    return true;
  } else 
    return false;
}

template <typename I, typename F>
__device__
inline bool static_magnetic_map(
    const I m2n0,
    const I nphi, // const I iphi, const I vphi,
    const I *m2tris,
    const F *m2invdet,
    const F *staticB,
    const bvh2d_node_t<I, F> *bvh,
    const F h, // delta,
    const I nsteps,
    F rz[2]) // input/output rz
{
  F v[2];
  for (int k = 0; k < nsteps; k ++) {
    if (!eval_static_b(m2n0, nphi, m2tris, m2invdet, staticB, bvh, rz, v))
      return false;
    const F k1[2] = {h * v[0], h * v[1]};

    const F rz2[2] = {rz[0] + k1[0]/2, rz[1] + k1[1]/2};
    if (!eval_static_b(m2n0, nphi, m2tris, m2invdet, staticB, bvh, rz2, v)) 
      return false;
    const F k2[2] = {h * v[0], h * v[1]};

    const F rz3[2] = {rz[0] + k2[0]/2, rz[1] + k2[1]/2};
    if (!eval_static_b(m2n0, nphi, m2tris, m2invdet, staticB, bvh, rz3, v)) 
      return false;
    const F k3[2] = {h * v[0], h * v[1]};
    
    const F rz4[2] = {rz[0] + k3[0], rz[1] + k3[1]};
    if (!eval_static_b(m2n0, nphi, m2tris, m2invdet, staticB, bvh, rz4, v)) 
      return false;
    const F k4[2] = {h * v[0], h * v[1]};
    
    for (int i = 0; i < 2; i ++) 
      rz[i] = rz[i] + (k1[i] + 2.0*(k2[i]+k3[i]) + k4[i])/6.0;
  }
  return true;
}

template <typename I, typename F>
__global__ void mx2_derive_interpolants(
    const I m2n0,
    const I nphi, const I iphi, const I vphi,
    const I *m2tris, const F *m2coords,
    const F *m2invdet,
    const F *staticB,
    const bvh2d_node_t<I, F> *bvh,
    const I p, // ith poloidal plane in vphi
    ftk::xgc_interpolant_t<I, F> *interpolants) // output
{
  int tid = getGlobalIdx_3D_1D();
  I i = tid;
  if (i >= m2n0) return;

  const I np = nphi * vphi;
  const I nip = nphi * iphi * vphi;
  const I steps_per_vp = 32; // TODO
  const F delta = 2 * M_PI / (nip * steps_per_vp);

  const I bsteps = p * steps_per_vp, 
          fsteps = (vphi - p) * steps_per_vp;

  F rz0[2] = { m2coords[i*2], m2coords[i*2+1] };
  ftk::xgc_interpolant_t<I, F> &interpolant = interpolants[i];

  { // backward
    F rz[2] = {rz0[0], rz0[1]};
    if ( static_magnetic_map(
        m2n0, nphi, // iphi, vphi, 
        m2tris, m2invdet, staticB, bvh, 
        -delta, bsteps, rz) ) {
      const I tid = bvh2_locate_point_tri(bvh, rz[0], rz[1], interpolant.mu0, m2invdet);
      if (tid >= 0) {
        m2_get_tri(tid, interpolant.tri0, m2tris);
      } else {
        interpolant.tri0[0] = -1;
      }
    }
  }
  
  { // forward 
    F rz[2] = {rz0[0], rz0[1]};
    if ( static_magnetic_map(
        m2n0, nphi, // iphi, vphi, 
        m2tris, m2invdet, staticB, bvh, 
        delta, fsteps, rz) ) {
      const I tid = bvh2_locate_point_tri(bvh, rz[0], rz[1], interpolant.mu1, m2invdet);
      if (tid >= 0) {
        m2_get_tri(tid, interpolant.tri1, m2tris);
      } else {
        interpolant.tri1[0] = -1;
      }
    }
  }
}

template <bool UseStaticB, typename I, typename F>
__device__
inline bool poincare_eval( // evaluate total B for computing poincare plot
    const I m2n0,
    // const I nphi, const I iphi, const I vphi,
    const I *m2tris,
    const F *m2invdet, 
    const F *staticB,
    /* const */ F **deltaB,
    const bvh2d_node_t<I, F> *bvh, 
    const F rz[2], // input rz
    const I p, // input poloidal plane
    F v[2]) // output normalized vector
{
  F mu[3]; // barycentric coordinates
  // const I tid = bvh2_locate_point_recursive(bvh, bvh, rz[0], rz[1], mu, m2invdet);
  const I tid = bvh2_locate_point_tri(bvh, rz[0], rz[1], mu, m2invdet);

  if (tid >= 0) {
    F B[3][3], b[3];
    I tri[3];
    m2_get_tri(tid, tri, m2tris);
    
    for (int i = 0; i < 3; i ++)
      for (int j = 0; j < 3; j ++) 
        if (UseStaticB)
          B[i][j] = staticB[tri[i]*3 + j]; 
        else 
          B[i][j] = staticB[tri[i]*3 + j] + deltaB[p][tri[i]*3 + j]; // + deltaB[p*m2n0*3 + tri[i]*3 + j];

    lerp_s2v3(B, mu, b);
    v[0] = rz[0] * b[0] / b[2];
    v[1] = rz[0] * b[1] / b[2];
    
    // printf("tid=%d, rz=%f, %f, v=%f, %f\n", tid, rz[0], rz[1], v[0], v[1]);
    return true;
  } else 
    return false;
}

template <bool UseStaticB, typename I, typename F>
__device__
inline bool poincare_integrate2pi_butcher_rk5(
    const I m2n0,
    const I nphi, const I iphi, const I vphi,
    const I *m2tris,
    const F *m2invdet, 
    const F *staticB,
    /* const */ F **deltaB,
    const bvh2d_node_t<I, F> *bvh,
    const int dir, // +1 or -1
    F rz[2]) // input/output rz
{
  const I np = nphi * vphi;
  const I nip = np * iphi;
  const F half_h = 2 * M_PI / nip * dir;
  const F quarter_h = half_h * 0.5;
  const F h = half_h * 2;

  F v[2]; 
  for (int ii = 0; ii < iphi; ii ++) { // for iphi
    for (int ip = 0; ip < np; ip += 4) {
      int p = dir == 1 ? ip : (np-ip);

      const F rz0[2] = {rz[0], rz[1]};
        
      if (!poincare_eval<UseStaticB, I, F>(m2n0, m2tris, m2invdet, staticB, deltaB, bvh, rz, p%np, v)) 
        return false;
      const F k1[2] = {v[0], v[1]};
      
      // printf("rz=%f, %f, v=%f, %f\n", rz[0], rz[1], v[0], v[1]);

      const F rz2[2] = {rz[0] + h*k1[0]/4, 
                        rz[1] + h*k1[1]/4};
      if (!poincare_eval<UseStaticB, I, F>(m2n0, m2tris, m2invdet, staticB, deltaB, bvh, rz2, p+dir, v))  // 1/4 delta
        return false;
      const F k2[2] = {v[0], v[1]};

      const F rz3[2] = {rz[0] + h*k1[0]/8 + h*k2[0]/8, 
                        rz[1] + h*k1[1]/8 + h*k2[1]/8};
      if (!poincare_eval<UseStaticB, I, F>(m2n0, m2tris, m2invdet, staticB, deltaB, bvh, rz3, p+dir, v)) // 1/4 delta
        return false;
      const F k3[2] = {v[0], v[1]};
      
      const F rz4[2] = {rz[0] - h*k2[0]/2 + h*k3[0], 
                        rz[1] - h*k2[1]/2 + h*k3[1]};
      if (!poincare_eval<UseStaticB, I, F>(m2n0, m2tris, m2invdet, staticB, deltaB, bvh, rz4, p+2*dir, v)) // 1/2 delta
        return false;
      const F k4[2] = {v[0], v[1]};
      
      const F rz5[2] = {rz[0] + h*k1[0]*3/16 + h*k4[0]*9/16, 
                        rz[1] + h*k1[1]*3/16 + h*k4[1]*9/16};
      if (!poincare_eval<UseStaticB, I, F>(m2n0, m2tris, m2invdet, staticB, deltaB, bvh, rz5, p+3*dir, v)) // 3/4 delta
        return false;
      const F k5[2] = {v[0], v[1]};
      
      const F rz6[2] = {rz[0] - h*k1[0]*3/7 + h*k2[0]*2/7 + h*k3[0]*12/7 - h*k4[0]*12/7 + h*k5[0]*8/7, 
                        rz[1] - h*k1[1]*3/7 + h*k2[1]*2/7 + h*k3[1]*12/7 - h*k4[1]*12/7 + h*k5[1]*8/7};
      if (!poincare_eval<UseStaticB, I, F>(m2n0, m2tris, m2invdet, staticB, deltaB, bvh, rz6, (p+4*dir)%np, v)) // 1 delta
        return false;
      const F k6[2] = {v[0], v[1]};

      for (int i = 0; i < 2; i ++) 
        rz[i] = rz[i] + h * (7*k1[i] + 32*k3[i] + 12*k4[i] + 32*k5[i] + 7*k6[i]) / 90;
    }
  }

  return true;
}

#if 0
template <bool UseStaticB, typename I, typename F>
__device__
inline bool poincare_integrate2pi(
    const I m2n0,
    const I nphi, const I iphi, const I vphi,
    const I *m2tris,
    const F *m2invdet, 
    const F *staticB,
    const F *deltaB,
    const bvh2d_node_t<I, F> *bvh, 
    F rz[2]) // input/output rz
{
  const I np = nphi * vphi;
  const F half_h = 2 * M_PI / np;
  const F h = half_h * 2;
  
  F v[2]; 
  for (int p = 0; p < np; p += 2) {
    const F rz0[2] = {rz[0], rz[1]};
      
    if (!poincare_eval<UseStaticB, I, F>(m2n0, nphi, iphi, vphi, m2tris, m2invdet, staticB, deltaB, bvh, rz, p, v)) 
      return false;
    const F k1[2] = {h * v[0], h * v[1]};
    
    // printf("rz=%f, %f, v=%f, %f\n", rz[0], rz[1], v[0], v[1]);

    const F rz2[2] = {rz[0] + k1[0]/2, rz[1] + k1[1]/2};
    if (!poincare_eval<UseStaticB, I, F>(m2n0, nphi, iphi, vphi, m2tris, m2invdet, staticB, deltaB, bvh, rz2, p+1, v)) 
      return false;
    const F k2[2] = {h * v[0], h * v[1]};

    const F rz3[2] = {rz[0] + k2[0]/2, rz[1] + k2[1]/2};
    if (!poincare_eval<UseStaticB, I, F>(m2n0, nphi, iphi, vphi, m2tris, m2invdet, staticB, deltaB, bvh, rz3, p+1, v)) 
      return false;
    const F k3[2] = {h * v[0], h * v[1]};
    
    const F rz4[2] = {rz[0] + k3[0], rz[1] + k3[1]};
    if (!poincare_eval<UseStaticB, I, F>(m2n0, nphi, iphi, vphi, m2tris, m2invdet, staticB, deltaB, bvh, rz4, (p+2)%np, v)) 
      return false;
    const F k4[2] = {h * v[0], h * v[1]};

    for (int i = 0; i < 2; i ++) 
      // rz[i] = rz[i] + (k1[i] + 2.0*(k2[i]+k3[i]) + h*v[i])/6.0;
      rz[i] = rz[i] + (k1[i] + 2.0*(k2[i]+k3[i]) + k4[i])/6.0;
  }

  return true;
}
#endif

template <typename I, typename F>
__global__ void poincare_compute_psin(
    const I m2n0,
    const I *m2tris,
    const F *m2invdet,
    const F *psin,
    const bvh2d_node_t<I, F> *bvh, 
    const I nseeds,
    const F *plot, 
    F *plot_psin)
{
  I i = getGlobalIdx_3D_1D();
  if (i >= nseeds) return;

  const F rz[2] = {plot[i*2], plot[i*2+1]};

  F mu[3];
  const I tid = bvh2_locate_point_tri(bvh, rz[0], rz[1], mu, m2invdet);

  if (tid >= 0) {
    F psins[3];
    I tri[3];
    m2_get_tri(tid, tri, m2tris);

    for (int j = 0; j < 3; j ++)
      psins[j] = psin[tri[j]];

    plot_psin[i] = lerp_s2(psins, mu);
  } else 
    plot_psin[i] = 1e38;
}

template <bool UseStaticB, typename I, typename F>
__global__ void poincare_integrate(
    const int k, // kth step
    const I m2n0, 
    const I nphi, const I iphi, const I vphi, 
    const I *m2tris,
    const F *m2invdet,
    const F *staticB,
    /* const */ F **deltaB, 
    const bvh2d_node_t<I, F> *bvh, 
    const I nseeds, 
    const I nsteps,
    const int dir,
    F *plot)
{
  int tid = getGlobalIdx_3D_1D();
  I i = tid;
  if (i >= nseeds) return;

  const auto offset = (k-1)*nseeds + i;
  F rz[2] = {plot[offset*2], plot[offset*2+1]};

  // bool succ = poincare_integrate2pi<UseStaticB, I, F>(
  bool succ = poincare_integrate2pi_butcher_rk5<UseStaticB, I, F>(
      m2n0, nphi, iphi, vphi,
      m2tris, m2invdet,
      staticB, deltaB, bvh, dir, rz);
  
  const auto offset1 = k*nseeds + i;
  if (succ) {
    plot[offset1*2] = rz[0];
    plot[offset1*2+1] = rz[1];
  } else {
    plot[offset1*2] = 1e38;
    plot[offset1*2+1] = 1e38;
    return;
  }
}

template <typename I, typename F>
__device__
inline bool check_simplex(
    I current_timestep, 
    I i,
    const F factor,
    const I nphi, const I iphi, const I vphi,
    const I m2n0, const I m2n1, const I m2n2,
    const F m2coords[], // m2 vertex coordinates
    const I m2edges[], // list of m2 edges
    const I m2tris[], // list of m2 triangles
    const F psin_[], // normalized psi
    const ftk::xgc_interpolant_t<I, F>* interpolants, // interpolants
    const F *const scalar[2], // current and next scalar
    const F *const vector[2], // current and next grad
    const F *const jacobian[2], // current and next jacobian
    cp_t & cp) // WIP: critical points
{
  // typedef ftk::fixed_point<> fp_t;

  const I np = nphi * iphi * vphi;
  const I m3n0 = m2n0 * np;

  I verts[3], t[3], p[3];
  F rzpt[3][4], f[3], v[3][2], j[3][2][2], psin[3];

  mx4_get_tri(i, verts, np, m2n0, m2n1, m2n2, m2edges, m2tris);

  for (int k = 0; k < 3; k ++) {
    t[k] = verts[k] / m3n0; // time in m4
    const I v3 = verts[k] % m3n0; // vert in m3
    const I v2 = v3 % m2n0; // vert in m2
    p[k] = v3 / m2n0; // poloidal plane
    psin[k] = psin_[v2];

    mx3_get_coords(v3, rzpt[k], m2n0, m2coords);
    rzpt[k][3] = t[k];
    // const int iv = (t[k] == current_timestep) ? 0 : 1;
    const int iv = t[k]; 
    mx3_interpolate<I, F>(
        v3, nphi, iphi, vphi, 
        m2n0, interpolants,
        scalar[iv], vector[iv], jacobian[iv], 
        f[k], v[k], j[k]);
  }

#if 0
  if (i > 100000000 && i < 100000100)
    printf("i=%d, verts=%d, %d, %d, t=%d, %d, %d, f=%f, %f, %f, p=%d, %d, %d, rzpt=%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
        i, verts[0], verts[1], verts[2], 
        t[0], t[1], t[2],
        f[0], f[1], f[2],
        p[0], p[1], p[2],
        rzpt[0][0], rzpt[0][1], rzpt[0][2], rzpt[0][3], 
        rzpt[1][0], rzpt[1][1], rzpt[1][2], rzpt[1][3],
        rzpt[2][0], rzpt[2][1], rzpt[2][2], rzpt[2][3]);
  return false;
#endif
  // printf("%f, %f, %f\n", f[0], f[1], f[2]);

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

  int64_t vf[3][2];
  for (int k = 0; k < 3; k ++) 
    for (int l = 0; l < 2; l ++) 
      vf[k][l] = factor * v[k][l];
 
  // if (v[0][0] > 0)
  //   printf("v=%f, %f, %f, %f, %f, %f\n", 
  //       v[0][0], v[0][1], v[1][0], v[1][1], v[2][0], v[2][1]);

  // if (vf[0][0] + vf[0][1] + vf[1][0] + vf[1][1] + vf[2][0] + vf[2][1] != 0)
  //   printf("vf=%lld, %lld; %lld, %lld; %lld, %lld\n", 
  //       vf[0][0], vf[0][1], vf[1][0], vf[1][1], vf[2][0], vf[2][1]);

  bool succ = ftk::robust_critical_point_in_simplex2(vf, verts);
  if (!succ) return false;

  F mu[3], x[4];
  bool succ2 = ftk::inverse_lerp_s2v2(v, mu);
  ftk::clamp_barycentric<3>(mu);

  ftk::lerp_s2v4(rzpt, mu, x);
  for (int k = 0; k < 3; k ++)
    cp.x[k] = x[k];
  cp.t = x[3];

  // cp.scalar[0] = f[0] * mu[0] + f[1] * mu[1] + f[2] * mu[2];
  cp.scalar[0] = ftk::lerp_s2( f, mu );
  cp.scalar[1] = ftk::lerp_s2( psin, mu );
  cp.tag = i;

  F h[2][2];
  ftk::lerp_s2m2x2(j, mu, h);
  // ftk::make_symmetric2x2(h);
  cp.type = ftk::critical_point_type_2d(h, true);
 
  // printf("tag=%d, verts=%d, %d, %d, x=%f, %f, %f, %f, scalar=%f\n", 
  //     i, verts[0], verts[1], verts[2], cp.x[0], cp.x[1], cp.x[2], cp.x[3], cp.scalar[0]);
  // printf("cp.x=%f, %f, %f, %f, scalar=%f, type=%d\n", 
  //     cp.x[0], cp.x[1], cp.x[2], cp.x[3], cp.scalar[0], 
  //     cp.type);

  return true;
}

template <typename I, typename F>
__global__
void sweep_simplices(
    int scope, 
    I current_timestep, 
    const F factor,
    const I nphi, const I iphi, const I vphi,
    const I m2n0, const I m2n1, const I m2n2,
    const F m2coords[], // m2 vertex coordinates
    const I m2edges[], // list of m2 edges
    const I m2tris[], // list of m2 triangles
    const F psin[],
    const ftk::xgc_interpolant_t<I, F> *interpolants, 
    const F* scalar0, // current scalar
    const F* scalar1, // next scalar
    const F* vector0, 
    const F* vector1, 
    const F* jacobian0, 
    const F* jacobian1,
    unsigned long long &ncps, 
    cp_t *cps)
{
  const I np = nphi * iphi * vphi;
  const I mx3n0 = m2n0 * np; 
  const I mx3n1 = (2 * m2n1 + m2n0) * np; 
  const I mx3n2 = (3 * m2n2 + 2 * m2n1) * np;
  const I mx4n2 = 3 * mx3n2 + 2 * mx3n1;

  int tid = getGlobalIdx_3D_1D();
  I i = tid;
  if (scope == scope_interval) {
    i += mx3n2;
    if (i >= mx4n2) return; // invalid element
  } else { // ordinal
    if (i >= mx3n2) return;
  }
 
  const F* const scalar[2] = {scalar0, scalar1};
  const F* const vector[2] = {vector0, vector1};
  const F* const jacobian[2] = {jacobian0, jacobian1};

  cp_t cp;
  bool succ = check_simplex<I, F>(
      current_timestep, 
      i, 
      factor,
      nphi, iphi, vphi, 
      m2n0, m2n1, m2n2, 
      m2coords, m2edges, m2tris,
      psin,
      interpolants, 
      scalar, vector, jacobian, 
      cp);
  
  if (succ) {
    unsigned long long idx = atomicAdd(&ncps, 1ul);
    cp.tag = i; // tid;
    cps[idx] = cp;
  }
}

void xft_create_poincare_ctx(ctx_t **c_, int nseeds, int nsteps, int device)
{
  xft_create_ctx(c_, device, std::ceil((nseeds * nsteps * sizeof(double) * 2.0) / (1024 * 1024))+16);
  ctx_t *c = *c_;
 
  c->nseeds = nseeds;
  c->nsteps = nsteps;
  c->poincare = true;
}

void xft_create_ctx(ctx_t **c_, int device, int device_buffer_size_in_mb)
{
  *c_ = (ctx_t*)malloc(sizeof(ctx_t));
  ctx_t *c = *c_;
  memset(c, 0, sizeof(ctx_t));

  c->device = device;
  cudaSetDevice(device);

  cudaMalloc((void**)&c->dncps, sizeof(unsigned long long));
  cudaMemset(c->dncps, 0, sizeof(unsigned long long));
  checkLastCudaError("[FTK-CUDA] cuda malloc");

  fprintf(stderr, "allocating %d MB\n", device_buffer_size_in_mb);
  c->bufsize = device_buffer_size_in_mb * size_t(1024 * 1024); 
  c->hcps = (cp_t*)malloc(c->bufsize);
  cudaMalloc((void**)&c->dcps, c->bufsize);
  checkLastCudaError("[FTK-CUDA] cuda malloc: creating buffer");

  c->d_psin = NULL;

  c->d_kernel_nodes = NULL;
  c->d_kernel_values = NULL;
  c->d_kernel_lengths = NULL;
  c->d_kernel_offsets = NULL;

  c->d_scalar_in = NULL;
  c->d_scalar[0] = NULL;
  c->d_scalar[1] = NULL;
  c->d_vector[0] = NULL;
  c->d_vector[1] = NULL;
  c->d_jacobian[0] = NULL;
  c->d_jacobian[1] = NULL;

  c->factor = 32768.0;
}

void xft_destroy_ctx(ctx_t **c_)
{
  ctx_t *c = *c_;

  if (c->d_m2coords != NULL) cudaFree(c->d_m2coords);
  if (c->d_m2edges != NULL) cudaFree(c->d_m2edges);
  if (c->d_m2tris != NULL) cudaFree(c->d_m2tris);
  if (c->d_psin != NULL) cudaFree(c->d_psin);
  if (c->d_m2invdet != NULL) cudaFree(c->d_m2invdet);

  if (c->d_vertex_triangles != NULL) cudaFree(c->d_vertex_triangles);

  if (c->d_bvh != NULL) cudaFree(c->d_bvh);

  if (c->d_interpolants != NULL) cudaFree(c->d_interpolants);

  if (c->d_kernel_nodes != NULL) cudaFree(c->d_kernel_nodes);
  if (c->d_kernel_values != NULL) cudaFree(c->d_kernel_values);
  if (c->d_kernel_lengths != NULL) cudaFree(c->d_kernel_lengths);
  if (c->d_kernel_offsets != NULL) cudaFree(c->d_kernel_offsets);

  if (c->d_scalar_in != NULL) cudaFree(c->d_scalar_in);
  if (c->d_scalar[0] != NULL) cudaFree(c->d_scalar[0]);
  if (c->d_scalar[1] != NULL) cudaFree(c->d_scalar[1]);
  if (c->d_vector[0] != NULL) cudaFree(c->d_vector[0]);
  if (c->d_vector[1] != NULL) cudaFree(c->d_vector[1]);
  if (c->d_jacobian[0] != NULL) cudaFree(c->d_jacobian[0]);
  if (c->d_jacobian[1] != NULL) cudaFree(c->d_jacobian[1]);
 
  if (c->d_apars != NULL) cudaFree(c->d_apars);
  if (c->d_apars_upsample != NULL) cudaFree(c->d_apars_upsample);
  if (c->d_bfield != NULL) cudaFree(c->d_bfield);
  if (c->d_bfield0 != NULL) cudaFree(c->d_bfield0);
  if (c->d_curl_bfield0 != NULL) cudaFree(c->d_curl_bfield0);
  // if (c->d_seeds != NULL) cudaFree(c->d_seeds);
  if (c->d_poincare_psin != NULL) cudaFree(c->d_poincare_psin);

  checkLastCudaError("[FTK-CUDA] cuda free");

  if (c->hcps != NULL) free(c->hcps);
  if (c->h_poincare_psin != NULL) free(c->h_poincare_psin);

  free(*c_);
  *c_ = NULL;
}

void xft_compute_poincare_psin(ctx_t *c)
{
  fprintf(stderr, "poincare_psin...\n");
  if (c->d_poincare_psin == NULL) {
    cudaMalloc((void**)&c->d_poincare_psin, sizeof(double) * c->nseeds * c->nsteps);
    checkLastCudaError("[FTK-CUDA] cuda compute poincare psin: malloc");
  }
  if (c->h_poincare_psin == NULL)
    c->h_poincare_psin = (double*)malloc(sizeof(double) * c->nseeds * c->nsteps);
  
  const int maxGridDim = 1024;
  const int blockSize = 256;
  const int nBlocks = idivup(c->nseeds * c->nsteps, blockSize);
  dim3 gridSize;
  if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else gridSize = dim3(nBlocks);

  poincare_compute_psin<<<gridSize, blockSize>>>(
      c->m2n0, c->d_m2tris, c->d_m2invdet, 
      c->d_psin, c->d_bvh, c->nseeds * c->nsteps, 
      (const double*)c->dcps, c->d_poincare_psin);
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] cuda compute poincare psin: kernel");

  cudaMemcpy(c->h_poincare_psin, c->d_poincare_psin, 
      sizeof(double) * c->nseeds * c->nsteps, 
      cudaMemcpyDeviceToHost);
  checkLastCudaError("[FTK-CUDA] cuda compute poincare psin: memcpy");
}

void xft_compute_poincare_plot(ctx_t *c, const double *seeds, bool use_static_b, int dir)
{
  // fprintf(stderr, "loading seeds..., %p, %p\n", c->dcps, seeds);
  cudaMemcpy(c->dcps, seeds, c->nseeds * sizeof(double) * 2, cudaMemcpyHostToDevice);
  // fprintf(stderr, "seeds loaded.\n");
  checkLastCudaError("[FTK-CUDA] xft_compute_poincare_plot: load seeds");
  
  const int maxGridDim = 1024;
  const int blockSize = 32;
  const int nBlocks = idivup(c->nseeds, blockSize);
  dim3 gridSize;

  if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else gridSize = dim3(nBlocks);

  for (int k = 1; k < c->nsteps; k ++) {
    if (k % 100 == 1) {
      cudaDeviceSynchronize();
      fprintf(stderr, "step k=%d, total steps %d\n", k, c->nsteps);
    }

    if (use_static_b)
      poincare_integrate<true, int, double><<<gridSize, blockSize>>>(
          k, c->m2n0, c->nphi, c->iphi, c->vphi, 
          c->d_m2tris, c->d_m2invdet, c->d_bfield, c->ddp_deltaB, 
          c->d_bvh, c->nseeds, c->nsteps, dir, (double*)c->dcps);
    else 
      poincare_integrate<false, int, double><<<gridSize, blockSize>>>(
          k, c->m2n0, c->nphi, c->iphi, c->vphi, 
          c->d_m2tris, c->d_m2invdet, c->d_bfield, c->ddp_deltaB, 
          c->d_bvh, c->nseeds, c->nsteps, dir, (double*)c->dcps);
  }
  checkLastCudaError("[FTK-CUDA] xft_compute_poincare_plot: exec");

  cudaMemcpy(c->hcps, c->dcps, c->nseeds * c->nsteps * sizeof(double) * 2, cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] xft_compute_poincare_plot: memcpy");
}

void xft_execute(ctx_t *c, int scope, int current_timestep)
{
  const int np = c->nphi * c->iphi * c->vphi;
  const int mx3n1 = (2 * c->m2n1 + c->m2n0) * np;
  const int mx3n2 = (3 * c->m2n2 + 2 * c->m2n1) * np;
  // const int mx4n2 = 3 * mx3n2 + 2 * mx3n1;
  const int mx4n2_ordinal  = mx3n2, 
            mx4n2_interval = 2 * mx3n2 + 2 * mx3n1;
  // fprintf(stderr, "executing timestep %d\n", current_timestep);

  size_t ntasks;
  if (scope == scope_ordinal) ntasks = mx4n2_ordinal;
  else ntasks = mx4n2_interval;
  
  fprintf(stderr, "ntasks=%zu\n", ntasks);
  
  const int maxGridDim = 1024;
  const int blockSize = 256;
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;

  if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else gridSize = dim3(nBlocks);

  sweep_simplices<int, double><<<gridSize, blockSize>>>(
      scope, current_timestep, 
      c->factor,
      c->nphi, c->iphi, c->vphi, 
      c->m2n0, c->m2n1, c->m2n2, 
      c->d_m2coords, c->d_m2edges, c->d_m2tris, 
      c->d_psin,
      c->d_interpolants, 
      c->d_scalar[0], c->d_scalar[1],
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
}

void xft_swap(ctx_t *c)
{
  // fprintf(stderr, "calling swap\n");
  std::swap(c->d_scalar[0], c->d_scalar[1]);
  std::swap(c->d_vector[0], c->d_vector[1]);
  std::swap(c->d_jacobian[0], c->d_jacobian[1]);
}

void xft_load_data(ctx_t *c, 
    const double *scalar, const double *vector, const double *jacobian)
{
  double *dd_scalar;
  if (c->d_scalar[0] == NULL) {
    cudaMalloc((void**)&c->d_scalar[0], sizeof(double) * size_t(c->m2n0) * size_t(c->nphi));
    checkLastCudaError("[FTK-CUDA] loading scalar field data, malloc 0");
    dd_scalar = c->d_scalar[0];
  } else if (c->d_scalar[1] == NULL) {
    cudaMalloc((void**)&c->d_scalar[1], sizeof(double) * size_t(c->m2n0) * size_t(c->nphi));
    checkLastCudaError("[FTK-CUDA] loading scalar field data, malloc 0.1");
    dd_scalar = c->d_scalar[1];
  } else {
    std::swap(c->d_scalar[0], c->d_scalar[1]);
    dd_scalar = c->d_scalar[1];
  }
  // fprintf(stderr, "dd=%p, d0=%p, d1=%p, src=%p\n", dd_scalar, c->d_scalar[0], c->d_scalar[1], scalar);
  cudaMemcpy(dd_scalar, scalar, sizeof(double) * size_t(c->m2n0 * c->nphi), 
      cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] loading scalar field data, memcpy 0");
 
  /// 
  double *dd_vector;
  if (c->d_vector[0] == NULL) {
    cudaMalloc((void**)&c->d_vector[0], sizeof(double) * size_t(c->m2n0) * size_t(c->nphi) * 2);
    dd_vector = c->d_vector[0];
  } else if (c->d_vector[1] == NULL) {
    cudaMalloc((void**)&c->d_vector[1], sizeof(double) * size_t(c->m2n0) * size_t(c->nphi)* 2);
    dd_vector = c->d_vector[1];
  } else {
    std::swap(c->d_vector[0], c->d_vector[1]);
    dd_vector = c->d_vector[1];
  }
  cudaMemcpy(dd_vector, vector, sizeof(double) * size_t(c->m2n0 * c->nphi * 2), 
      cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] loading vector field data");

  /// 
  double *dd_jacobian;
  if (c->d_jacobian[0] == NULL) {
    cudaMalloc((void**)&c->d_jacobian[0], sizeof(double) * size_t(c->m2n0) * size_t(c->nphi) * 4);
    dd_jacobian = c->d_jacobian[0];
  } else if (c->d_jacobian[1] == NULL) {
    cudaMalloc((void**)&c->d_jacobian[1], sizeof(double) * size_t(c->m2n0) * size_t(c->nphi) * 4);
    dd_jacobian = c->d_jacobian[1];
  } else {
    std::swap(c->d_jacobian[0], c->d_jacobian[1]);
    dd_jacobian = c->d_jacobian[1];
  }
  cudaMemcpy(dd_jacobian, jacobian, sizeof(double) * size_t(c->m2n0 * c->nphi) * 4, 
      cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] loading jacobian field data");
}

template <typename I, typename F>
__global__
void smooth_scalar_vector_jacobian(
    const I nphi,
    const I m2n0,
    const F *m2coords,
    const F sigma,
    const size_t *lengths,
    const size_t *offsets,
    const I *nodes,
    const F *values,
    const F *scalar_in,
    F *scalar_out,
    F *vector_out, 
    F *jacobian_out)
{
  int idx = getGlobalIdx_3D_1D();
  if (idx >= m2n0 * nphi) return; // out of bounds

  const F sigma2 = sigma * sigma, 
          sigma4 = sigma2 * sigma2;

  const I i = idx % m2n0, 
          p = idx / m2n0;

  // scalar_out[idx] = 0; // assuming out values are zero'ed
  for (int j = 0; j < lengths[i]; j ++) {
    const int k = nodes[offsets[i] + j];
    const F w = values[offsets[i] + j]; // weight
    const F d[2] = {m2coords[k*2] - m2coords[i*2], m2coords[k*2+1] - m2coords[i*2+1]};
  
    const F f = scalar_in[ k + p * m2n0 ];
    scalar_out[idx] += f * w;
    
    // if (idx == 10000)
    // if (scalar_in[idx] > 0.01)
    //   printf("w[%d]=%f, scalar_in=%f\n", k, w, scalar_in[idx]);

    vector_out[idx*2]   += -f * w * d[0] / sigma2;
    vector_out[idx*2+1] += -f * w * d[1] / sigma2;

    jacobian_out[idx*4]   += (d[0]*d[0] / sigma2 - 1) / sigma2 * f * w;
    jacobian_out[idx*4+1] += d[0]*d[1] / sigma4 * f * w;
    jacobian_out[idx*4+2] += d[0]*d[1] / sigma4 * f * w;
    jacobian_out[idx*4+3] += (d[1]*d[1] / sigma2 - 1) / sigma2 * f * w;
  }
    
  // if (scalar_out[idx] > 0.01)
  //   printf("%f, %f, %f\n", scalar_out[idx], vector_out[idx*2], vector_out[idx*2+1]);
}

void xft_smooth_scalar_vector_jacobian(ctx_t *c, 
    const double *d_scalar_in, 
    double *d_scalar_out, 
    double *d_vector_out, 
    double *d_jacobian_out)
{
  cudaMemset(d_scalar_out, 0, size_t(c->m2n0 * c->nphi) * sizeof(double));
  cudaMemset(d_vector_out, 0, size_t(2 * c->m2n0 * c->nphi) * sizeof(double));
  cudaMemset(d_jacobian_out, 0, size_t(4 * c->m2n0 * c->nphi) * sizeof(double));
  
  size_t ntasks = c->nphi * c->m2n0;
  const int maxGridDim = 1024;
  const int blockSize = 256;
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;
  if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else gridSize = dim3(nBlocks);

  smooth_scalar_vector_jacobian<int, double><<<gridSize, blockSize>>>(
      c->nphi, 
      c->m2n0, 
      c->d_m2coords,
      c->sigma,
      c->d_kernel_lengths,
      c->d_kernel_offsets,
      c->d_kernel_nodes,
      c->d_kernel_values,
      d_scalar_in, 
      d_scalar_out,
      d_vector_out,
      d_jacobian_out);
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] smoothing scalar vector jacobian");

  // update scaling factor; CPU implementation for now..
  double *h_vector = (double*)malloc(size_t(c->m2n0 * c->nphi) * sizeof(double));
  cudaMemcpy(h_vector, d_vector_out, size_t(c->m2n0 * c->nphi) * sizeof(double), 
      cudaMemcpyDeviceToHost);
  double maxabs = 0.0;
  for (int i = 0; i < c->m2n0 * c->nphi; i ++)
    maxabs = std::max(maxabs, std::abs(h_vector[i]));
  free(h_vector);
  
  double factor = std::exp2(-std::ceil(std::log2(maxabs)) + 20); // 20 bits
  c->factor = std::max(c->factor, factor);

  std::cerr << "maxabs: " << maxabs << ", factor: " << c->factor << std::endl;
}

void xft_load_smoothing_kernel(ctx_t *c, double sigma, const std::vector<std::vector<std::tuple<int, double>>>& kernels)
{
  c->sigma = sigma;

  // fprintf(stderr, "loading smoothing kernels to GPU...\n");
  std::vector<size_t> lengths(kernels.size()), offsets(kernels.size());
  std::vector<int> nodes;
  std::vector<double> values;

  size_t acc = 0;
  for (size_t i = 0; i < kernels.size(); i ++) { // nodes
    const std::vector<std::tuple<int, double>>& kernel = kernels[i];
    lengths[i] = kernel.size();
    offsets[i] = acc;
    acc += kernel.size();

    for (size_t j = 0; j < kernel.size(); j ++) {
      nodes.push_back(std::get<0>(kernel[j]));
      values.push_back(std::get<1>(kernel[j]));
    }
  }
  
  cudaMalloc((void**)&c->d_kernel_nodes, nodes.size() * sizeof(int));
  cudaMemcpy(c->d_kernel_nodes, nodes.data(), nodes.size() * sizeof(int), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_kernel_values, values.size() * sizeof(double));
  cudaMemcpy(c->d_kernel_values, values.data(), values.size() * sizeof(double), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_kernel_lengths, lengths.size() * sizeof(size_t));
  cudaMemcpy(c->d_kernel_lengths, lengths.data(), lengths.size() * sizeof(size_t), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_kernel_offsets, offsets.size() * sizeof(size_t));
  cudaMemcpy(c->d_kernel_offsets, offsets.data(), offsets.size() * sizeof(size_t), cudaMemcpyHostToDevice);

  checkLastCudaError("[FTK-CUDA] loading smoothing kernel");
  // fprintf(stderr, "smoothing kernels loaded to GPU.\n");
}

void xft_load_psin(ctx_t *c, const double *psin)
{
  cudaMalloc((void**)&c->d_psin, size_t(c->m2n0) * sizeof(double));
  cudaMemcpy(c->d_psin, psin, size_t(c->m2n0) * sizeof(double), cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] loading psin");
}

void xft_load_magnetic_field(ctx_t *c, const double *bfield, const double *bfield0, const double *curl_bfield0)
{
  cudaMalloc((void**)&c->d_bfield, size_t(c->m2n0) * sizeof(double) * 3);
  cudaMemcpy(c->d_bfield, bfield, size_t(c->m2n0) * sizeof(double) * 3, cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_bfield0, size_t(c->m2n0) * sizeof(double) * 3);
  cudaMemcpy(c->d_bfield0, bfield0, size_t(c->m2n0) * sizeof(double) * 3, cudaMemcpyHostToDevice);
  
  cudaMalloc((void**)&c->d_curl_bfield0, size_t(c->m2n0) * sizeof(double) * 3);
  cudaMemcpy(c->d_curl_bfield0, curl_bfield0, size_t(c->m2n0) * sizeof(double) * 3, cudaMemcpyHostToDevice);

  checkLastCudaError("[FTK-CUDA] loading xgc magnetic field");
}

void xft_load_apars(ctx_t *c, const double *apars)
{
  const int maxGridDim = 1024;
  const int blockSize = 256;
  
  if (c->d_apars == NULL)
    cudaMalloc((void**)&c->d_apars, size_t(c->m2n0 * c->nphi) * sizeof(double));
  cudaMemcpy(c->d_apars, apars, size_t(c->m2n0 * c->nphi) * sizeof(double), cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] xft_load_apars: copy apars");

  cudaDeviceSynchronize();
  fprintf(stderr, "apars copied\n");

  // TODO:
  // 1. upsample apars
  if (c->d_apars_upsample == NULL) {
    cudaMalloc((void**)&c->d_apars_upsample, size_t(c->m2n0) * c->nphi * c->vphi * sizeof(double));
  }
  {
    const int nBlocks = idivup(c->m2n0 * c->nphi * c->vphi, blockSize);
    dim3 gridSize;
    if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
    else gridSize = dim3(nBlocks);

    mx3_upsample_scalar<<<gridSize, blockSize>>>(
        c->nphi, c->iphi, c->vphi, c->m2n0,
        c->d_interpolants,
        c->d_apars, c->d_apars_upsample);

    cudaFree(c->d_apars);
    c->d_apars = NULL;
  }
  checkLastCudaError("[FTK-CUDA] xft_load_apars: upsample apars");
  
  cudaDeviceSynchronize();
  fprintf(stderr, "apars upsampled\n");
 
  if (c->retrieve_apars_upsample) {
    // if (c->h_apars_upsample)
    c->h_apars_upsample = (double*)malloc(c->m2n0 * c->nphi * c->vphi * sizeof(double));

    cudaMemcpy(c->h_apars_upsample, c->d_apars_upsample,
        c->m2n0 * c->nphi * c->vphi * sizeof(double), 
        cudaMemcpyDeviceToHost);

    fprintf(stderr, "apars_upsample retrieved.\n");
  }

#if 0
  // 2. derive gradient of upsampled apars
  if (c->d_gradAs == NULL) {
    // cudaMalloc((void**)&c->d_gradAs_cw, 2 * size_t(c->m2n2) * c->nphi * c->vphi * sizeof(double));
    cudaMalloc((void**)&c->d_gradAs_cw, 2 * size_t(c->m2n2) * sizeof(double));
    cudaDeviceSynchronize();
    checkLastCudaError("[FTK-CUDA] xft_load_apars: malloc gradAs_cw");
    cudaMalloc((void**)&c->d_gradAs, 2 * size_t(c->m2n0) * c->nphi * c->vphi * sizeof(double));
    cudaDeviceSynchronize();
    checkLastCudaError("[FTK-CUDA] xft_load_apars: malloc gradAs");
  }
  for (int i = 0; i < c->nphi * c->vphi; i ++) {
    // fprintf(stderr, "deriving gradAs for plane %d\n", i);
    {
      const int nBlocks = idivup(c->m2n2, blockSize);
      dim3 gridSize;
      if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
      else gridSize = dim3(nBlocks);

      m2_cellwise_scalar_gradient<<<gridSize, blockSize>>>(
          c->m2n2, c->d_m2tris, c->d_m2coords,
          c->d_apars_upsample + i * c->m2n0, 
          // c->d_gradAs_cw + i * c->m2n2 * 2);
          c->d_gradAs_cw);
    
      cudaDeviceSynchronize();
      // fprintf(stderr, "cw\n");
      checkLastCudaError("[FTK-CUDA] xft_load_apars: cw");
    }

    {
      const int nBlocks = idivup(c->m2n0, blockSize);
      dim3 gridSize;
      if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
      else gridSize = dim3(nBlocks);

      // fprintf(stderr, "mt=%d\n", c->max_vertex_triangles);
      m2_cellwise2vertexwise_scalar_gradient<<<gridSize, blockSize>>>(
          c->m2n0, c->max_vertex_triangles, c->d_vertex_triangles, 
          c->d_gradAs_cw, c->d_gradAs + i * c->m2n0 * 2);
      
      cudaDeviceSynchronize();
      // fprintf(stderr, "vw\n");
      checkLastCudaError("[FTK-CUDA] xft_load_apars: vw");
    }
  }
    
  cudaFree(c->d_gradAs_cw);
  c->d_gradAs_cw = NULL;
  
  fprintf(stderr, "grad_apars derived\n");
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] xft_load_apars: grad of upsampled apars");
#endif

  // 3. derive deltaB
#if 0
  if (c->d_deltaB == NULL) {
    fprintf(stderr, "allocating deltaB, %llu\n", 3ll * c->m2n0 * c->nphi * c->vphi * sizeof(double));
    cudaMalloc((void**)&c->d_deltaB, 3ll * c->m2n0 * c->nphi * c->vphi * sizeof(double));
    cudaDeviceSynchronize();
    checkLastCudaError("[FTK-CUDA] xft_load_apars: derive deltaB: malloc");
  }
#endif
  {
    const int nBlocks = idivup(std::max(c->m2n0, c->m2n2), blockSize);
    dim3 gridSize;
    if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
    else gridSize = dim3(nBlocks);

    // partial delta B
    double *d_partial_deltaB, *h_full_deltaB = (double*)malloc(3 * c->nphi * c->vphi * c->m2n0 * sizeof(double));
    cudaMalloc(&d_partial_deltaB, 3 * c->m2n0 * sizeof(double));

    double *d_partial_gradAs_cw; //  = d_partial_deltaB; // reusing the memory for gradAs_cw
    cudaMalloc(&d_partial_gradAs_cw, 2 * c->m2n2 * sizeof(double));
    double *d_partial_gradAs;
    cudaMalloc(&d_partial_gradAs, 2 * c->m2n0 * sizeof(double));

    for (int p = 0; p < c->nphi * c->vphi; p ++) { // the full deltaB is too large; compute each individually and copy back to host before copy to device again
      // fprintf(stderr, "computing deltaB for plane %d\n", p);
      // cudaMalloc((void**)&c->hdp_deltaB[p], 3 * c->m2n0 * sizeof(double));

      // cellwise gradient
      m2_cellwise_scalar_gradient<<<gridSize, blockSize>>>(
          c->m2n2, c->d_m2tris, c->d_m2coords, 
          c->d_apars_upsample + p * c->m2n0, 
          d_partial_gradAs_cw);

      // vertexwise gradient
      m2_cellwise2vertexwise_scalar_gradient<<<gridSize, blockSize>>>(
          c->m2n0, c->max_vertex_triangles, c->d_vertex_triangles, 
          d_partial_gradAs_cw, d_partial_gradAs);

      // fprintf(stderr, "deriving deltaB for plane %d..\n", p);
      mx3_derive_deltaB<<<gridSize, blockSize>>>(
          c->nphi, c->iphi, c->vphi, c->m2n0,
          c->d_apars_upsample,
          d_partial_gradAs, // c->d_gradAs,
          c->d_bfield,
          c->d_bfield0,
          c->d_curl_bfield0,
          p,
          d_partial_deltaB); // c->hdp_deltaB[p]); // c->d_deltaB);

      cudaMemcpy(h_full_deltaB + p * 3 * c->m2n0, 
          d_partial_deltaB, 
          3 * c->m2n0 * sizeof(double), 
          cudaMemcpyDeviceToHost);
    }
    cudaDeviceSynchronize();
    checkLastCudaError("[FTK-CUDA] xft_load_apars: derive deltaB: kernel");

    cudaFree(d_partial_gradAs_cw);
    cudaFree(d_partial_gradAs);
    cudaFree(d_partial_deltaB);
    cudaFree(c->d_gradAs); c->d_gradAs = NULL;
    cudaFree(c->d_apars_upsample); c->d_apars_upsample = NULL;
    
    cudaDeviceSynchronize();
    checkLastCudaError("[FTK-CUDA] xft_load_apars: free up partial arrays");
    
    c->hdp_deltaB = (double**)malloc(c->nphi * c->vphi * sizeof(double*));
    for (int p = 0; p < c->nphi * c->vphi; p ++) {
      // fprintf(stderr, "copying back deltaB to plane %d\n", p);
      cudaMalloc((void**)&c->hdp_deltaB[p], 3 * c->m2n0 * sizeof(double));
      checkLastCudaError("[FTK-CUDA] xft_load_apars: allocating partial deltaB");
      cudaMemcpy(c->hdp_deltaB[p], 
          h_full_deltaB + p * 3 * c->m2n0,
          3 * c->m2n0 * sizeof(double),
          cudaMemcpyHostToDevice);
      cudaDeviceSynchronize();
      checkLastCudaError("[FTK-CUDA] xft_load_apars: copy deltaB back to gpu");
    }
    free(h_full_deltaB);

    cudaMalloc((void**)&c->ddp_deltaB, c->nphi * c->vphi * sizeof(double*));
    cudaMemcpy(c->ddp_deltaB, c->hdp_deltaB, 
        c->nphi * c->vphi * sizeof(double*), 
        cudaMemcpyHostToDevice);
  }
  
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] xft_load_apars: derive deltaB");
  
  fprintf(stderr, "deltaB derived\n");
}

#if 0
void xft_retrieve_deltaB(ctx_t *c)
{
  if (c->h_deltaB == NULL)
    c->h_deltaB = (double*)malloc(c->m2n0 * c->nphi * c->vphi * sizeof(double) * 3);

  cudaMemcpy(c->h_deltaB, c->d_deltaB, 
      c->m2n0 * c->nphi * c->vphi * sizeof(double) * 3, 
      cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();

  checkLastCudaError("[FTK-CUDA] xft_retrieve_deltaB");
  fprintf(stderr, "deltaB retrieved\n");
}
#endif

void xft_load_vertex_triangles(ctx_t *c, const std::vector<std::set<int>> &vertex_triangles)
{
  size_t mt = 0;
  for (const auto &tris : vertex_triangles)
    mt = std::max(mt, tris.size());

  std::vector<int> h_vertex_triangles(mt * c->m2n0, -1);
  for (int i = 0; i < c->m2n0; i ++) {
    int j = 0;
    for (const auto tri : vertex_triangles[i]) {
      h_vertex_triangles[i*mt + j] = tri;
      j ++;
    }
  }

  // fprintf(stderr, "#vertex=%zu, mt=%zu\n", vertex_triangles.size(), mt);

  if (c->d_vertex_triangles != NULL) cudaFree(c->d_vertex_triangles);
  cudaMalloc((void**)&c->d_vertex_triangles, size_t(c->m2n0 * mt) * sizeof(int));
  cudaMemcpy(c->d_vertex_triangles, &h_vertex_triangles[0], size_t(c->m2n0 * mt) * sizeof(int), cudaMemcpyHostToDevice);
  c->max_vertex_triangles = mt;
  checkLastCudaError("[FTK-CUDA] xft_load_vertex_triangles");
}

void xft_load_bvh(ctx_t *c, const std::vector<bvh2d_node_t<>>& bvh)
{
  cudaMalloc((void**)&c->d_bvh, bvh.size() * sizeof(bvh2d_node_t<>));
  cudaMemcpy(c->d_bvh, &bvh[0], bvh.size() * sizeof(bvh2d_node_t<>), cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] xft_load_bvh");
  fprintf(stderr, "bvh loaded.\n");
  // cudaDeviceSynchronize();
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

  cudaMalloc((void**)&c->d_m2coords, size_t(m2n0) * sizeof(double) * 2);
  checkLastCudaError("[FTK-CUDA] loading xgc mesh, malloc 0");
  cudaMemcpy(c->d_m2coords, m2coords, size_t(m2n0) * sizeof(double) * 2, cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] loading xgc mesh, memcpy 0");

  cudaMalloc((void**)&c->d_m2edges, size_t(m2n1) * sizeof(int) * 2);
  cudaMemcpy(c->d_m2edges, m2edges, size_t(m2n1) * sizeof(int) * 2, cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_m2tris, size_t(m2n2) * sizeof(int) * 3);
  cudaMemcpy(c->d_m2tris, m2tris, size_t(m2n2) * sizeof(int) * 3, cudaMemcpyHostToDevice);
 
  if (c->poincare) { // compute invdet for point locator
    cudaMalloc((void**)&c->d_m2invdet, size_t(m2n2) * sizeof(double));
  
    const int maxGridDim = 1024;
    const int blockSize = 256;
    const int nBlocks = idivup(c->m2n2, blockSize);
    dim3 gridSize;
    if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
    else gridSize = dim3(nBlocks);

    fprintf(stderr, "computing m2invdet..\n");
    m2_compute_invdet<<<gridSize, blockSize>>>(
        c->m2n2, c->d_m2tris, c->d_m2coords,
        c->d_m2invdet);
  }

  checkLastCudaError("[FTK-CUDA] loading xgc mesh");
}

void xft_derive_interpolants(ctx_t *c)
{
  if (c->d_interpolants == NULL)
    cudaMalloc((void**)&c->d_interpolants,
        size_t(c->m2n0) * sizeof(ftk::xgc_interpolant_t<>) * (c->vphi-1));
      
  const int maxGridDim = 1024;
  const int blockSize = 256;
  const int nBlocks = idivup(c->m2n0, blockSize);
  dim3 gridSize;
  if (nBlocks >= maxGridDim) gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else gridSize = dim3(nBlocks);
  
  cudaDeviceSynchronize();

  for (size_t p = 1; p < c->vphi; p ++) {
  fprintf(stderr, "deriving interpolants for plane %d...\n", p);
    mx2_derive_interpolants<int, double><<<gridSize, blockSize>>>(
        c->m2n0, c->nphi, c->iphi, c->vphi, 
        c->d_m2tris, c->d_m2coords, c->d_m2invdet, 
        c->d_bfield, c->d_bvh, p,
        c->d_interpolants + (p-1) * size_t(c->m2n0));
    cudaDeviceSynchronize();
  }

  fprintf(stderr, "interpolants derived.\n");
  checkLastCudaError("[FTK-CUDA] deriving interpolants");
}

void xft_load_interpolants(ctx_t *c, const std::vector<std::vector<ftk::xgc_interpolant_t<>>> &interpolants)
{
  fprintf(stderr, "loading interpolants, %zu, %zu\n", c->vphi, interpolants.size());
  assert(c->vphi == interpolants.size());

  cudaMalloc((void**)&c->d_interpolants, 
      size_t(c->m2n0) * sizeof(ftk::xgc_interpolant_t<>) * c->vphi);
  checkLastCudaError("[FTK-CUDA] loading xgc interpolants, malloc 0");

  for (size_t i = 1; i < interpolants.size(); i ++) {
    cudaMemcpy(c->d_interpolants + (i-1) * size_t(c->m2n0), // * sizeof(ftk::xgc_interpolant_t<>), 
        interpolants[i].data(), size_t(c->m2n0) * sizeof(ftk::xgc_interpolant_t<>), cudaMemcpyHostToDevice);
    checkLastCudaError("[FTK-CUDA] loading xgc interpolants, memcpy");
  }
  
  checkLastCudaError("[FTK-CUDA] loading xgc interpolants");
  
  // fprintf(stderr, "interpolants loaded.\n");
  // cudaDeviceSynchronize();
}

void xft_load_scalar_data(ctx_t *c, const double *scalar)
{
  // fprintf(stderr, "loading and smoothing scalar data w/ gpu...\n");
  if (c->d_scalar_in == NULL)
    cudaMalloc((void**)&c->d_scalar_in, sizeof(double) * size_t(c->m2n0 * c->nphi));
  cudaMemcpy(c->d_scalar_in, scalar, sizeof(double) * size_t(c->m2n0 * c->nphi), 
      cudaMemcpyHostToDevice);

  double *dd_scalar, *dd_vector, *dd_jacobian;
  if (c->d_scalar[0] == NULL) {
    // fprintf(stderr, "init slot 0\n");
    cudaMalloc((void**)&c->d_scalar[0], sizeof(double) * size_t(c->m2n0 * c->nphi));
    cudaMalloc((void**)&c->d_vector[0], sizeof(double) * size_t(c->m2n0 * c->nphi * 2));
    cudaMalloc((void**)&c->d_jacobian[0], sizeof(double) * size_t(c->m2n0 * c->nphi * 4));
    dd_scalar = c->d_scalar[0];
    dd_vector = c->d_vector[0];
    dd_jacobian = c->d_jacobian[0];
  } else if (c->d_scalar[1] == NULL) {
    // fprintf(stderr, "init slot 1\n");
    cudaMalloc((void**)&c->d_scalar[1], sizeof(double) * size_t(c->m2n0 * c->nphi));
    cudaMalloc((void**)&c->d_vector[1], sizeof(double) * size_t(c->m2n0 * c->nphi * 2));
    cudaMalloc((void**)&c->d_jacobian[1], sizeof(double) * size_t(c->m2n0 * c->nphi * 4));
    dd_scalar = c->d_scalar[1];
    dd_vector = c->d_vector[1];
    dd_jacobian = c->d_jacobian[1];
  } else {
    // fprintf(stderr, "swapping 0 and 1\n");
    std::swap(c->d_scalar[0], c->d_scalar[1]);
    std::swap(c->d_vector[0], c->d_vector[1]);
    std::swap(c->d_jacobian[0], c->d_jacobian[1]);
    dd_scalar = c->d_scalar[1];
    dd_vector = c->d_vector[1];
    dd_jacobian = c->d_jacobian[1];
  }

  xft_smooth_scalar_vector_jacobian(c, 
      c->d_scalar_in, dd_scalar, dd_vector, dd_jacobian);
  
  // fprintf(stderr, "scalar smoothed and loaded to gpu\n");
}
