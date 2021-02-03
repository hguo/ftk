#include <CL/sycl.hpp>
#include <vector>
#include <ftk/config.hh>
#include <ftk/mesh/lattice.hh>
#include "common.cuh"

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
  std::vector<cp_t> results;
  const size_t ntasks = core.n() * ntypes_3_2<scope>();

  cl::sycl::buffer<double> dVc(Vc, 2 * sizeof(double) * ext.n());

  return results;
}

std::vector<cp_t> extract_cp2dt_sycl(
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
