#ifndef _FTK_LEVELSET_TRACKER
#define _FTK_LEVELSET_TRACKER

#include <ftk/filters/filter.hh>

namespace ftk {

struct levelset_tracker : public filter
{
  levelset_tracker() {}
  virtual ~levelset_tracker();
};

}

#endif
