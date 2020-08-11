#ifndef _FTK_FILTER_HH
#define _FTK_FILTER_HH

#include <ftk/ftk_config.hh>
#include <ftk/object.hh>
#include <ftk/external/cxxopts.hpp>
#include <thread>
#include <mutex>
#include <thread>
#include <cassert>

namespace ftk {

enum { 
  FTK_XL_NONE = 0,
  FTK_XL_OPENMP = 1,
  FTK_XL_SYCL = 2,
  FTK_XL_TBB = 3,
  FTK_XL_CUDA = 4,
  FTK_XL_KOKKOS_CUDA = 5
};

struct filter : public object {
  filter() {
    nthreads = default_nthreads();
  }

  virtual void update() = 0;
  virtual void reset() {};

  void use_accelerator(int i) {
    xl = i;
    if (xl == FTK_XL_OPENMP || xl == FTK_XL_SYCL || xl == FTK_XL_TBB || xl == FTK_XL_KOKKOS_CUDA) {
      warn("Accelerator not available.  Using FTK_XL_NONE.");
      xl = FTK_XL_NONE;
    }
  }

  int default_nthreads() const {
    if (comm.size() > 1) return 1; // use 1 thread per proc for mpi runs
    else return std::thread::hardware_concurrency(); 
  }

  void set_communicator(const diy::mpi::communicator comm_) {comm = comm_;}
  
  void set_number_of_threads(int n) {nthreads = n;}
  int get_number_of_threads() const {return nthreads;}

  void set_root_proc(int p) {root_proc = p;}
  int get_root_proc() const {return root_proc;}
  bool is_root_proc() const {return root_proc == comm.rank();}

  void set_device_id(int d) {set_device_ids(std::vector<int>({d}));}
  void set_device_ids(const std::vector<int>& ids) {device_ids = ids;}
  const std::vector<int>& get_device_ids() const {return device_ids;}
  int get_number_devices() const {return device_ids.size();}

protected:
  diy::mpi::communicator comm;

  int xl = FTK_XL_NONE;
  int nthreads = 1;
  int root_proc = 0;
  std::vector<int> device_ids;

  std::mutex mutex;
};

}

#endif
