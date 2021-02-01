#include <nvfunctional>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/clamp.hh>
#include <ftk/numeric/symmetric_matrix.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/filters/xgc_blob_filament_tracker.cuh>
#include "common.cuh"
#include "mx4.cuh"

// what are needed in device memory:
// - field-following interpolants for virtual poloidal planes
// - triangles, edges, and vertex coordinates in 2D mesh
// - current and next timestep of scalar, vector, and jacobian fields

using namespace ftk;

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
    const F psin_[], // normalized psi
    const ftk::xgc_interpolant_t<I, F>* interpolants, // interpolants
    const F *const scalar[2], // current and next scalar
    const F *const vector[2], // current and next grad
    const F *const jacobian[2], // current and next jacobian
    cp_t & cp) // WIP: critical points
{
  // typedef ftk::fixed_point<> fp_t;
  const long long factor = 1 << 15; // WIP

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

  long long vf[3][2];
  for (int k = 0; k < 3; k ++) 
    for (int l = 0; l < 2; l ++) 
      vf[k][l] = factor * v[k][l];
  
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
    const I nphi, const I iphi, const I vphi,
    const I m2n0, const I m2n1, const I m2n2,
    const F m2coords[], // m2 vertex coordinates
    const I m2edges[], // list of m2 edges
    const I m2tris[], // list of m2 triangles
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
      nphi, iphi, vphi, 
      m2n0, m2n1, m2n2, 
      m2coords, m2edges, m2tris,
      interpolants, 
      scalar, vector, jacobian, 
      cp);
  
  if (succ) {
    unsigned long long idx = atomicAdd(&ncps, 1ul);
    cp.tag = i; // tid;
    cps[idx] = cp;
  }
}

void xft_create_ctx(ctx_t **c_)
{
  *c_ = (ctx_t*)malloc(sizeof(ctx_t));
  ctx_t *c = *c_;
  
  cudaMalloc((void**)&c->dncps, sizeof(unsigned long long));
  cudaMemset(c->dncps, 0, sizeof(unsigned long long));
  checkLastCudaError("[FTK-CUDA] cuda malloc");

  c->bufsize = 512 * 1024 * 1024; 
  c->hcps = (cp_t*)malloc(c->bufsize);
  cudaMalloc((void**)&c->dcps, c->bufsize);
  checkLastCudaError("[FTK-CUDA] cuda malloc");

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
}

void xft_destroy_ctx(ctx_t **c_)
{
  ctx_t *c = *c_;

  if (c->d_m2coords != NULL) cudaFree(c->d_m2coords);
  if (c->d_m2edges != NULL) cudaFree(c->d_m2edges);
  if (c->d_m2tris != NULL) cudaFree(c->d_m2tris);

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
  
  checkLastCudaError("[FTK-CUDA] cuda free");

  free(*c_);
  *c_ = NULL;
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
      c->nphi, c->iphi, c->vphi, 
      c->m2n0, c->m2n1, c->m2n2, 
      c->d_m2coords, c->d_m2edges, c->d_m2tris, 
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
    cudaMalloc((void**)&c->d_scalar[0], sizeof(double) * c->m2n0 * c->nphi);
    checkLastCudaError("[FTK-CUDA] loading scalar field data, malloc 0");
    dd_scalar = c->d_scalar[0];
  } else if (c->d_scalar[1] == NULL) {
    cudaMalloc((void**)&c->d_scalar[1], sizeof(double) * c->m2n0 * c->nphi);
    checkLastCudaError("[FTK-CUDA] loading scalar field data, malloc 0.1");
    dd_scalar = c->d_scalar[1];
  } else {
    std::swap(c->d_scalar[0], c->d_scalar[1]);
    dd_scalar = c->d_scalar[1];
  }
  // fprintf(stderr, "dd=%p, d0=%p, d1=%p, src=%p\n", dd_scalar, c->d_scalar[0], c->d_scalar[1], scalar);
  cudaMemcpy(dd_scalar, scalar, sizeof(double) * c->m2n0 * c->nphi, 
      cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] loading scalar field data, memcpy 0");
 
  /// 
  double *dd_vector;
  if (c->d_vector[0] == NULL) {
    cudaMalloc((void**)&c->d_vector[0], sizeof(double) * c->m2n0 * c->nphi * 2);
    dd_vector = c->d_vector[0];
  } else if (c->d_vector[1] == NULL) {
    cudaMalloc((void**)&c->d_vector[1], sizeof(double) * c->m2n0 * c->nphi * 2);
    dd_vector = c->d_vector[1];
  } else {
    std::swap(c->d_vector[0], c->d_vector[1]);
    dd_vector = c->d_vector[1];
  }
  cudaMemcpy(dd_vector, vector, sizeof(double) * c->m2n0 * c->nphi * 2, 
      cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] loading vector field data");

  /// 
  double *dd_jacobian;
  if (c->d_jacobian[0] == NULL) {
    cudaMalloc((void**)&c->d_jacobian[0], sizeof(double) * c->m2n0 * c->nphi * 4);
    dd_jacobian = c->d_jacobian[0];
  } else if (c->d_jacobian[1] == NULL) {
    cudaMalloc((void**)&c->d_jacobian[1], sizeof(double) * c->m2n0 * c->nphi * 4);
    dd_jacobian = c->d_jacobian[1];
  } else {
    std::swap(c->d_jacobian[0], c->d_jacobian[1]);
    dd_jacobian = c->d_jacobian[1];
  }
  cudaMemcpy(dd_jacobian, jacobian, sizeof(double) * c->m2n0 * c->nphi * 4, 
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
  cudaMemset(d_scalar_out, 0, c->m2n0 * c->nphi * sizeof(double));
  cudaMemset(d_vector_out, 0, 2 * c->m2n0 * c->nphi * sizeof(double));
  cudaMemset(d_jacobian_out, 0, 4 * c->m2n0 * c->nphi * sizeof(double));
  
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
  cudaMalloc((void**)&c->d_psin, m2n0 * sizeof(double));
  cudaMemcpy(c->d_psin, psin, m2n0 * sizeof(double), cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] loading psin");
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

  cudaMalloc((void**)&c->d_m2coords, m2n0 * sizeof(double) * 2);
  checkLastCudaError("[FTK-CUDA] loading xgc mesh, malloc 0");
  cudaMemcpy(c->d_m2coords, m2coords, m2n0 * sizeof(double) * 2, cudaMemcpyHostToDevice);
  checkLastCudaError("[FTK-CUDA] loading xgc mesh, memcpy 0");

  cudaMalloc((void**)&c->d_m2edges, m2n1 * sizeof(int) * 2);
  cudaMemcpy(c->d_m2edges, m2edges, m2n1 * sizeof(int) * 2, cudaMemcpyHostToDevice);

  cudaMalloc((void**)&c->d_m2tris, m2n2 * sizeof(int) * 3);
  cudaMemcpy(c->d_m2tris, m2tris, m2n2 * sizeof(int) * 3, cudaMemcpyHostToDevice);
  
  checkLastCudaError("[FTK-CUDA] loading xgc mesh");
}

void xft_load_interpolants(ctx_t *c, const std::vector<std::vector<ftk::xgc_interpolant_t<>>> &interpolants)
{
  assert(c->vphi == interpolants.size());

  cudaMalloc((void**)&c->d_interpolants, 
      c->m2n0 * sizeof(ftk::xgc_interpolant_t<>) * c->vphi);
  checkLastCudaError("[FTK-CUDA] loading xgc interpolants, malloc 0");

  for (int i = 1; i < interpolants.size(); i ++) {
    cudaMemcpy(c->d_interpolants + (i-1) * c->m2n0, // * sizeof(ftk::xgc_interpolant_t<>), 
        interpolants[i].data(), c->m2n0 * sizeof(ftk::xgc_interpolant_t<>), cudaMemcpyHostToDevice);
    checkLastCudaError("[FTK-CUDA] loading xgc interpolants, memcpy");
  }
  
  checkLastCudaError("[FTK-CUDA] loading xgc interpolants");
}

void xft_load_scalar_data(ctx_t *c, const double *scalar)
{
  // fprintf(stderr, "loading and smoothing scalar data w/ gpu...\n");
  if (c->d_scalar_in == NULL)
    cudaMalloc((void**)&c->d_scalar_in, sizeof(double) * c->m2n0 * c->nphi);
  cudaMemcpy(c->d_scalar_in, scalar, sizeof(double) * c->m2n0 * c->nphi, 
      cudaMemcpyHostToDevice);

  double *dd_scalar, *dd_vector, *dd_jacobian;
  if (c->d_scalar[0] == NULL) {
    // fprintf(stderr, "init slot 0\n");
    cudaMalloc((void**)&c->d_scalar[0], sizeof(double) * c->m2n0 * c->nphi);
    cudaMalloc((void**)&c->d_vector[0], sizeof(double) * c->m2n0 * c->nphi * 2);
    cudaMalloc((void**)&c->d_jacobian[0], sizeof(double) * c->m2n0 * c->nphi * 4);
    dd_scalar = c->d_scalar[0];
    dd_vector = c->d_vector[0];
    dd_jacobian = c->d_jacobian[0];
  } else if (c->d_scalar[1] == NULL) {
    // fprintf(stderr, "init slot 1\n");
    cudaMalloc((void**)&c->d_scalar[1], sizeof(double) * c->m2n0 * c->nphi);
    cudaMalloc((void**)&c->d_vector[1], sizeof(double) * c->m2n0 * c->nphi * 2);
    cudaMalloc((void**)&c->d_jacobian[1], sizeof(double) * c->m2n0 * c->nphi * 4);
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
