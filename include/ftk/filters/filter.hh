#ifndef _FTK_FILTER_HH
#define _FTK_FILTER_HH

#include <ftk/ftk_config.hh>
#include <ftk/external/diy/mpi.hpp>
#include <mutex>

namespace ftk {

struct filter {
  filter() {
    if (comm.size() > 1) nthreads = 1;
    else nthreads = std::thread::hardware_concurrency(); // TODO
  }

  diy::mpi::communicator comm;

  int nthreads = 1;
  std::mutex mutex;

  virtual void update() = 0;
};

}

#endif
