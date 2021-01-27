#ifndef _FTK_FILTER_HH
#define _FTK_FILTER_HH

#include <ftk/config.hh>
#include <ftk/object.hh>
#include <ftk/external/cxxopts.hpp>
#include <thread>
#include <mutex>
#include <thread>
#include <cassert>

namespace ftk {

struct filter : public object {
  filter(diy::mpi::communicator comm) : object(comm) {
    nthreads = default_nthreads();
  }

  virtual void update() = 0;
  virtual void reset() {};

  void use_accelerator(const std::string& acc);
  void use_accelerator(int i) {
    xl = i;
#if 0
    if (xl == FTK_XL_OPENMP || xl == FTK_XL_SYCL || xl == FTK_XL_KOKKOS_CUDA) {
      warn("Accelerator not available.  Using FTK_XL_NONE.");
      xl = FTK_XL_NONE;
    }
#endif
  }

  int default_nthreads() const {
    if (comm.size() > 1) return 1; // use 1 thread per proc for mpi runs
    else return std::thread::hardware_concurrency(); 
  }
  
  void set_number_of_threads(int n) {nthreads = n;}
  int get_number_of_threads() const {return nthreads;}

  void set_number_of_blocks(int n) {nblocks = n; fprintf(stderr, "setting nb=%d\n", n);}
  int get_number_of_blocks() const {return nblocks;}

  void set_device_id(int d) {set_device_ids(std::vector<int>({d}));}
  void set_device_ids(const std::vector<int>& ids) {device_ids = ids;}
  const std::vector<int>& get_device_ids() const {return device_ids;}
  int get_number_devices() const {return device_ids.size();}

protected:
  int xl = FTK_XL_NONE;
  int nthreads = 1, nblocks = 0;
  bool enable_set_affinity = false; // true;

  std::vector<int> device_ids;

  std::mutex mutex;
};

////
inline void filter::use_accelerator(const std::string& acc)
{
  if (acc == "cuda") use_accelerator(FTK_XL_CUDA);
  else if (acc == "openmp") use_accelerator(FTK_XL_OPENMP);
  else if (acc == "tbb") use_accelerator(FTK_XL_TBB);
  else if (acc == "sycl") use_accelerator(FTK_XL_SYCL);
  else use_accelerator(FTK_XL_NONE);
}

}

#endif
