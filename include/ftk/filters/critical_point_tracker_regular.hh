#ifndef _FTK_CRITICAL_POINT_TRACKER_REGULAR_HH
#define _FTK_CRITICAL_POINT_TRACKER_REGULAR_HH

#include <ftk/filters/critical_point_tracker.hh>

namespace ftk {

struct critical_point_tracker_regular : public critical_point_tracker {
  critical_point_tracker_regular() {}
  critical_point_tracker_regular(int argc, char **argv) : critical_point_tracker(argc, argv) {}
};

}

#endif
