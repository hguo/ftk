#ifndef _FTK_FEATURE_SURFACE_HH
#define _FTK_FEATURE_SURFACE_HH

#include <ftk/filters/feature_curve.hh>

namespace ftk {

struct feature_surface_t {
  std::vector<feature_point_t> pts;
  std::vector<std::array<int, 3>> conn;
};

}

#endif
