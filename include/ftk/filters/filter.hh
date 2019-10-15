#ifndef _FTK_FILTER_HH
#define _FTK_FILTER_HH

#include <ftk/ftk_config.hh>
#include <mutex>

namespace ftk {

struct filter {
  int nthreads = 1;
  std::mutex mutex;
};

}

#endif
