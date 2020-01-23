#ifndef _FTK_CRITICAL_POINT_TEST_HH
#define _FTK_CRITICAL_POINT_TEST_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/sign_det.hh>

// reference:
// Edelsbrunner and Mucke, Simulation of simplicity: A technique to cope with degenerate cases in geometric algorithms.

namespace ftk {

template <typename FixedPointType=long long, typename WeightType=int>
__device__ __host__
inline bool robust_critical_point_in_simplex2(const FixedPointType V[3][2], const WeightType indices[3])
{
  const FixedPointType zero[2] = {0};
  return robust_point_in_simplex2(V, indices, zero, WeightType(-1)); // -1 is the index of zero point 
}

template <typename FixedPointType=long long, typename WeightType=int>
__device__ __host__
inline bool robust_critical_point_in_simplex3(const FixedPointType V[4][3], const WeightType indices[4])
{
  const FixedPointType zero[3] = {0};
  return robust_point_in_simplex3(V, indices, zero, WeightType(-1));
}

} // namespace ftk

#endif
