#ifndef _FTK_FILTER_HH
#define _FTK_FILTER_HH

#include <ftk/config.hh>
#include <ftk/object.hh>
#include <ftk/external/cxxopts.hpp>
#include <ndarray/util.hh>
#include <thread>
#include <mutex>
#include <thread>
#include <cassert>

namespace ftk {

/**
 * @brief Base class for all FTK filters
 *
 * The filter class provides the foundation for all data processing and analysis
 * operations in FTK. It handles parallelization (threading and MPI), GPU acceleration,
 * and provides a common interface for data transformation pipelines.
 *
 * Filters can be configured to use different execution backends (pthread, OpenMP, TBB)
 * and hardware accelerators (CUDA, SYCL). Thread and device management is handled
 * automatically based on the runtime environment.
 */
struct filter : public object {
  filter(diy::mpi::communicator comm) : object(comm) {
    nthreads = default_nthreads();
  }

  /**
   * @brief Process data through the filter
   *
   * This pure virtual function must be implemented by derived classes to define
   * the core data processing logic of the filter.
   */
  virtual void update() = 0;

  /**
   * @brief Reset the filter to its initial state
   *
   * Clears any cached data or internal state. Default implementation does nothing.
   */
  virtual void reset() {};

  /**
   * @brief Set the threading backend by name
   * @param backend Backend name: "openmp", "tbb", or "pthread" (default)
   */
  void use_thread_backend(const std::string& backend);

  /**
   * @brief Set the threading backend by ID
   * @param i Backend ID (FTK_THREAD_PTHREAD, FTK_THREAD_OPENMP, or FTK_THREAD_TBB)
   */
  void use_thread_backend(int i) { thread_backend = i; }

  /**
   * @brief Set the hardware accelerator by name
   * @param acc Accelerator name: "cuda", "sycl", or "none" (default)
   */
  void use_accelerator(const std::string& acc);

  /**
   * @brief Set the hardware accelerator by ID
   * @param i Accelerator ID (FTK_XL_NONE, FTK_XL_CUDA, or FTK_XL_SYCL)
   */
  void use_accelerator(int i) {
    xl = i;
#if 0
    if (xl == FTK_THREAD_OPENMP || xl == FTK_XL_SYCL || xl == FTK_XL_KOKKOS_CUDA) {
      warn("Accelerator not available.  Using FTK_XL_NONE.");
      xl = FTK_XL_NONE;
    }
#endif
  }

  /**
   * @brief Get the default number of threads for this environment
   * @return 1 if MPI is used with multiple processes, otherwise hardware_concurrency()
   */
  int default_nthreads() const {
    if (comm.size() > 1) return 1; // use 1 thread per proc for mpi runs
    else return std::thread::hardware_concurrency();
  }

  /**
   * @brief Set the number of threads for parallel operations
   * @param n Number of threads to use
   */
  void set_number_of_threads(int n) {nthreads = n;}

  /**
   * @brief Get the configured number of threads
   * @return Number of threads
   */
  int get_number_of_threads() const {return nthreads;}

  /**
   * @brief Set the number of blocks for domain decomposition
   * @param n Number of blocks
   */
  void set_number_of_blocks(int n) {nblocks = n; fprintf(stderr, "setting nb=%d\n", n);}

  /**
   * @brief Get the number of blocks
   * @return Number of blocks
   */
  int get_number_of_blocks() const {return nblocks;}

  /**
   * @brief Set a single GPU device ID
   * @param d Device ID
   */
  void set_device_id(int d) {set_device_ids(std::vector<int>({d}));}

  /**
   * @brief Set multiple GPU device IDs
   * @param ids Vector of device IDs to use
   */
  void set_device_ids(const std::vector<int>& ids) {device_ids = ids;}

  /**
   * @brief Set GPU device IDs from comma-separated string
   * @param ids Comma-separated device IDs (e.g., "0,1,2")
   */
  void set_device_ids(const std::string &ids);

  /**
   * @brief Get the list of GPU device IDs
   * @return Vector of device IDs
   */
  const std::vector<int>& get_device_ids() const {return device_ids;}

  /**
   * @brief Get the number of GPU devices configured
   * @return Number of devices
   */
  int get_number_devices() const {return device_ids.size();}

  /**
   * @brief Set the GPU device buffer size
   * @param mb Buffer size in megabytes
   */
  void set_device_buffer_size(int mb) { device_buffer_size_in_mb = mb; }

protected:
  int xl = FTK_XL_NONE, thread_backend = FTK_THREAD_PTHREAD;
  int nthreads = 1, nblocks = 0;
  bool enable_set_affinity = false; // true;

  std::vector<int> device_ids;
  int device_buffer_size_in_mb = 512;

  std::mutex mutex;
};

////
inline void filter::use_thread_backend(const std::string& str)
{
  if (str == "openmp") use_thread_backend( FTK_THREAD_OPENMP );
  else if (str == "tbb") use_thread_backend( FTK_THREAD_TBB );
  else use_thread_backend( FTK_THREAD_PTHREAD );
}

inline void filter::use_accelerator(const std::string& acc)
{
  if (acc == "cuda") use_accelerator(FTK_XL_CUDA);
  else if (acc == "sycl") use_accelerator(FTK_XL_SYCL);
  else use_accelerator(FTK_XL_NONE);
}

inline void filter::set_device_ids(const std::string& ids)
{
  if (ids.empty()) return;

  auto strs = split(ids, ",");
  std::vector<int> myids;
  for (auto str : strs)
    myids.push_back( std::stoi(str) );
  // for (auto i : myids) 
  //   fprintf(stderr, "using device %d\n", i);
  set_device_ids(myids);
}

}

#endif
