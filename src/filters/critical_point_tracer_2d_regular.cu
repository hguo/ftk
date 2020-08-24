#include <nvfunctional>
#include <cstdio>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/symmetric_matrix.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/mesh/lattice.hh>
#include <ftk/filters/critical_point.hh>
#include "common.cuh"

//// 
template <int scope>
__device__ __host__
bool check_simplex_cp2t(
    int current_timestep,
    const lattice3_t& domain,
    const lattice3_t& core, 
    const lattice2_t& ext, 
    const element32_t& e, 
    const double *V[2], // current and next timesteps
    const double *gradV[2], // jacobians
    const double *scalar[2], // scalars
    bool use_explicit_coords,
    const double *coords, // coordinates of vertices
    cp_t &cp)
{
  typedef ftk::fixed_point<> fp_t;

  // const int last_timestep = current_timestep - 1;
  // if (scope == scope_interval && e.corner[2] != current_timestep) // last_timestep)
  if (e.corner[2] != current_timestep) // last_timestep)
    return false;

  int vertices[3][3], indices[3];
  size_t local_indices[3];
  for (int i = 0; i < 3; i ++) {
    for (int j = 0; j < 3; j ++) {
      vertices[i][j] = e.corner[j] 
        + unit_simplex_offset_3_2<scope>(e.type, i, j);
      if (vertices[i][j] < domain.st[j] || 
           vertices[i][j] > domain.st[j] + domain.sz[j] - 1)
      // if (vertices[i][j] < core.st[j] || 
      //     vertices[i][j] > core.st[j] + core.sz[j] - 1)
        return false;
    }
    indices[i] = domain.to_index(vertices[i]);
    local_indices[i] = ext.to_index(vertices[i]);
  }

  double v[3][2];
  fp_t vf[3][2];
  for (int i = 0; i < 3; i ++) {
    // size_t k = ext.to_index(vertices[i]);
    const size_t k = local_indices[i];
    for (int j = 0; j < 2; j ++) {
      v[i][j] = V[unit_simplex_offset_3_2<scope>(e.type, i, 2/*time dimension id*/)][k*2+j];
      vf[i][j] = v[i][j];
    }
  }
  
  bool succ = robust_critical_point_in_simplex2(vf, indices);
  if (succ) {
    // inverse interpolation
    double mu[3];
    ftk::inverse_lerp_s2v2(v, mu, 0.0);
    // linear jacobian interpolation
    if (gradV[1]) { // have given jacobian
      double Js[3][2][2], J[2][2];
      for (int i = 0; i < 3; i ++) {
        // size_t ii = ext.to_index(vertices[i]);
        const size_t ii = local_indices[i];
        const int t = unit_simplex_offset_3_2<scope>(e.type, i, 2);
        for (int j = 0; j < 2; j ++) 
          for (int k = 0; k < 2; k ++)
            Js[i][j][k] = gradV[t][ii*4 + j*2 + k];
      }
      ftk::lerp_s2m2x2(Js, mu, J);
      ftk::make_symmetric2x2(J);
      cp.type = ftk::critical_point_type_2d(J, true/*symmetric*/);
      // if (cp.type != 0x100) return false;
    }

    // scalar interpolation
    if (scalar[0]) { // have given scalar
      double values[3];
      for (int i = 0; i < 3; i ++) {
        // const size_t ii = ext.to_index(vertices[i]);
        const size_t ii = local_indices[i];
        const int t = unit_simplex_offset_3_2<scope>(e.type, i, 2);
        values[i] = scalar[t][ii];
      }
      cp.scalar[0] = ftk::lerp_s2(values, mu);
      // if (abs(cp.scalar) < 0.02) return false; // threshold
    }
    
    double X[3][3], x[3];
#if 0 // TODO: use explicit coords
    if (use_explicit_coords) {
      for (int i = 0; i < 3; i ++) {
        for (int j = 0; j < 2; j ++) {
          X[i][j] = coords[j + local_indices[i]*2];
        }
        X[i][2] = vertices[i][2]; // unit_simplex_offset_3_2<scope>(e.type, i, 2);
      }
      ftk::lerp_s2v3(X, mu, cp.rx);
    }
#endif
    
    // implicit coordinates
    for (int i = 0; i < 3; i ++)
      for (int j = 0; j < 3; j ++)
        X[i][j] = vertices[i][j];
    ftk::lerp_s2v3(X, mu, x);
    cp.x[0] = x[0];
    cp.x[1] = x[1];
    cp.t = x[2];
    
    return true;
  } else 
    return false;
}

template <int scope>
__global__
void sweep_simplices(
    int current_timestep,
    const lattice3_t domain,
    const lattice3_t core,
    const lattice2_t ext, // array dimensions
    const double *Vc, // current timestep
    const double *Vn, // next timestep
    const double *Jc, 
    const double *Jn,
    const double *Sc, 
    const double *Sn,
    bool use_explicit_coords,
    const double *coords, // coordinates of vertices
    unsigned long long &ncps, cp_t *cps)
{
  const double *V[2] = {Vc, Vn};
  const double *J[2] = {Jc, Jn};
  const double *S[2] = {Sc, Sn};
  
  int tid = getGlobalIdx_3D_1D();
  const element32_t e = element32_from_index<scope>(core, tid);

  cp_t cp;
  bool succ = check_simplex_cp2t<scope>(
      current_timestep, 
      domain, core, ext, e, V, J, S, 
      use_explicit_coords, coords,
      cp);

  if (succ) {
    unsigned long long i = atomicAdd(&ncps, 1ul);
    cp.tag = tid;
    cps[i] = cp;
  }
}

template<int scope>
static std::vector<cp_t> extract_cp2dt(
    int current_timestep,
    const lattice3_t& domain,
    const lattice3_t& core, 
    const lattice2_t& ext, 
    const double *Vc, // 3D array: 2*W*H
    const double *Vn, 
    const double *Jc,
    const double *Jn,
    const double *Sc,
    const double *Sn,
    bool use_explicit_coords,
    const double *coords)
{
  const size_t ntasks = core.n() * ntypes_3_2<scope>();
  const int maxGridDim = 1024;
  const int blockSize = 256;
  const int nBlocks = idivup(ntasks, blockSize);
  dim3 gridSize;
  
  if (nBlocks >= maxGridDim) 
    gridSize = dim3(idivup(nBlocks, maxGridDim), maxGridDim);
  else 
    gridSize = dim3(nBlocks);

  double *dcoords = NULL;
  if (coords) {
    cudaMalloc((void**)&dcoords, 2 * sizeof(double) * ext.n());
    cudaMemcpy(dcoords, coords, 2 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }

  double *dVc, *dVn = NULL;
  if (Vc) {
    cudaMalloc((void**)&dVc, 2 * sizeof(double) * ext.n());
    cudaMemcpy(dVc, Vc, 2 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }
  if (Vn) {
    cudaMalloc((void**)&dVn, 2 * sizeof(double) * ext.n());
    cudaMemcpy(dVn, Vn, 2 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }

  double *dJc = NULL, *dJn = NULL;
  if (Jc) {
    cudaMalloc((void**)&dJc, 4 * sizeof(double) * ext.n());
    cudaMemcpy(dJc, Jc, 4 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }
  if (Jn) {
    cudaMalloc((void**)&dJn, 4 * sizeof(double) * ext.n());
    cudaMemcpy(dJn, Jn, 4 * sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }

  double *dSc = NULL, *dSn = NULL;
  if (Sc) {
    cudaMalloc((void**)&dSc, sizeof(double) * ext.n());
    cudaMemcpy(dSc, Sc, sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }
  if (Sn) {
    cudaMalloc((void**)&dSn, sizeof(double) * ext.n());
    cudaMemcpy(dSn, Sn, sizeof(double) * ext.n(), cudaMemcpyHostToDevice);
  }

  unsigned long long *dncps; // number of cps
  cudaMalloc((void**)&dncps, sizeof(unsigned long long));
  cudaMemset(dncps, 0, sizeof(unsigned long long));

  cp_t *dcps;
  // cudaMalloc((void**)&dcps, sizeof(cp_t) * ext.n() * 2);
  cudaMalloc((void**)&dcps, 1024*1024*64); // sizeof(cp_t) * ext.n() * 2);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMalloc/cudaMemcpy");

  fprintf(stderr, "calling kernel func...\n");
  sweep_simplices<scope><<<gridSize, blockSize>>>(
      current_timestep, 
      domain, core, ext, dVc, dVn, dJc, dJn, dSc, dSn,
      use_explicit_coords, dcoords, 
      *dncps, dcps);
  cudaDeviceSynchronize();
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices");

  unsigned long long ncps;
  cudaMemcpy(&ncps, dncps, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMemcpy, ncps");
  fprintf(stderr, "ncps=%llu\n", ncps);

  std::vector<cp_t> cps(ncps);
  cudaMemcpy(cps.data(), dcps, sizeof(cp_t) * ncps, cudaMemcpyDeviceToHost);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaMemcpy, dcps");
 
  if (dcoords) cudaFree(dcoords);
  if (dVc) cudaFree(dVc);
  if (dVn) cudaFree(dVn);
  if (dJc) cudaFree(dJc);
  if (dJn) cudaFree(dJn);
  if (dSc) cudaFree(dSc);
  if (dSn) cudaFree(dSn);
  cudaFree(dncps);
  cudaFree(dcps);
  checkLastCudaError("[FTK-CUDA] error: sweep_simplices: cudaFree");
 
  cudaDeviceSynchronize();
  fprintf(stderr, "exitting gpu kernel, ncps=%llu\n", ncps);

  return cps;
}

std::vector<cp_t>
extract_cp2dt_cuda(
    int scope, 
    int current_timestep,
    const ftk::lattice& domain,
    const ftk::lattice& core, 
    const ftk::lattice& ext, 
    const double *Vc, 
    const double *Vn, 
    const double *Jc, 
    const double *Jn, 
    const double *Sc,
    const double *Sn, 
    bool use_explicit_coords,
    const double *coords)
{
  lattice3_t D(domain);
  lattice3_t C(core);
  lattice2_t E(ext);

  // std::cerr << "domain=" << domain 
  //   << ", core=" << core << ", current_timestep=" 
  //   << current_timestep << std::endl;

  if (scope == scope_interval) 
    return extract_cp2dt<scope_interval>(current_timestep, 
        D, C, E, Vc, Vn, Jc, Jn, Sc, Sn, 
        use_explicit_coords, coords);
  if (scope == scope_ordinal) 
    return extract_cp2dt<scope_ordinal>(current_timestep, 
        D, C, E, Vc, Vn, Jc, Jn, Sc, Sn,
        use_explicit_coords, coords);
  else // scope == 2
    return extract_cp2dt<scope_all>(current_timestep, 
        D, C, E, Vc, Vn, Jc, Jn, Sc, Sn,
        use_explicit_coords, coords);
}
