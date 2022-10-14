#ifndef _FTK_FEATURE_POINT_LITE_HH
#define _FTK_FEATURE_POINT_LITE_HH

#include <ftk/config.hh>

namespace ftk {

struct feature_point_lite_t {
  double x[3] = {0}; // spatial coordinates 
  double t = 0.0; // time
  // double cond = 0.0; // condition number
  double scalar[FTK_CP_MAX_NUM_VARS] = {0};
  unsigned int type = 0;
  unsigned long long tag = 0;
};

}

#endif
