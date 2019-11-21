#ifndef _FTK_FILTER_HH
#define _FTK_FILTER_HH

#include <ftk/ftk_config.hh>
#include <ftk/external/diy/mpi.hpp>
#include <mutex>

namespace ftk {

enum { 
  FTK_XL_NONE = 0,
  FTK_XL_OPENMP,
  FTK_XL_SYCL,
  FTK_XL_TBB,
  FTK_XL_CUDA,
  FTK_XL_KOKKOS_CUDA
};

struct filter {
  filter() {
    if (comm.size() > 1) nthreads = 1; // use 1 thread per proc for mpi runs
    else nthreads = std::thread::hardware_concurrency(); 
  }

  virtual void update() = 0;

  void use_accelerator(int i) {xl = i;}

protected:
  diy::mpi::communicator comm;

  int xl = FTK_XL_NONE;
  int nthreads = 1;
  std::mutex mutex;
};

}

#endif
