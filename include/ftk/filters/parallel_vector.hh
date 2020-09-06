#ifndef _FTK_PARALLEL_VECTOR_T_HH
#define _FTK_PARALLEL_VECTOR_T_HH

#include <ftk/ftk_config.hh>
#include <array>
#include <ftk/external/json.hh>

namespace ftk {

using nlohmann::json;

struct parallel_vector_point_t {
  std::array<double, 3> x = {0};
  double t = 0.0;
  double lambda = 0.0;
  std::array<double, 3> v = {0}, w = {0}; // v = lambda * w
  double cond = 0.0;
  int timestep;
  bool ordinal = false;
  unsigned long long tag = 0;
};

}

#endif
