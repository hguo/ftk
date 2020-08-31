#ifndef _FTK_CRITICAL_POINT_LITE_HH
#define _FTK_CRITICAL_POINT_LITE_HH

#include <ftk/ftk_config.hh>

namespace ftk {

struct critical_point_lite_t {
  double x[3] = {0}; // spatial coordinates 
  double t = 0.0; // time
  double scalar[FTK_CP_MAX_NUM_VARS] = {0};
  unsigned int type = 0;
  unsigned long long tag = 0;
};

}

#endif
