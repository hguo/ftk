#ifndef _FTK_OBJECT_HH
#define _FTK_OBJECT_HH

#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <set>
#include <map>
#include <functional>
#include <algorithm>
#include <ftk/config.hh>
#include <ftk/error.hh>
#include <ftk/external/diy/mpi.hpp>
#include <unistd.h>
#include <sched.h>
#include <sys/resource.h>
#include <sys/syscall.h>

#if FTK_HAVE_TBB
#include <tbb/tbb.h>
#endif

namespace ftk {

/**
 * @brief Threading backend identifiers
 */
enum { // thread backend
  FTK_THREAD_NONE = 0,     ///< No threading (serial execution)
  FTK_THREAD_PTHREAD = 0,  ///< POSIX threads (default)
  FTK_THREAD_OPENMP = 1,   ///< OpenMP threading
  FTK_THREAD_TBB = 6       ///< Intel TBB threading
};

/**
 * @brief Hardware accelerator identifiers
 */
enum {
  FTK_XL_NONE = 0,         ///< No accelerator (CPU only)
  FTK_XL_SYCL = 2,         ///< SYCL acceleration
  FTK_XL_CUDA = 4,         ///< CUDA acceleration
  FTK_XL_KOKKOS_CUDA = 5   ///< Kokkos CUDA acceleration
};

/**
 * @brief Base class for all FTK objects
 *
 * This class provides common functionality for all FTK classes:
 * - MPI communicator management for distributed processing
 * - Root process designation for I/O operations
 * - Parallel loop execution with multiple threading backends
 * - CPU affinity control for performance tuning
 *
 * The object class serves as the foundation of FTK's execution model,
 * enabling both shared-memory and distributed-memory parallelism.
 */
struct object {
  object() {}
  object(diy::mpi::communicator c) {comm = c;}

  /**
   * @brief Set the MPI communicator
   * @param comm_ MPI communicator to use
   */
  void set_communicator(const diy::mpi::communicator comm_) {comm = comm_;}

  /**
   * @brief Set the root process for I/O operations
   * @param p Root process rank
   */
  void set_root_proc(int p) {root_proc = p;}

  /**
   * @brief Get the root process rank
   * @return Root process rank
   */
  int get_root_proc() const {return root_proc;}

  /**
   * @brief Check if this is the root process
   * @return True if this process is the root process
   */
  bool is_root_proc() const {return root_proc == comm.rank();}

  /**
   * @brief Set CPU affinity for the current thread
   * @param cpu CPU core ID
   */
  static void set_affinity(int cpu) {
#if !defined(_MSC_VER) && !defined(__APPLE__)
    cpu_set_t cpu_set;
    CPU_ZERO(&cpu_set);
    CPU_SET(cpu, &cpu_set);

    pthread_t thread = pthread_self();
    pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpu_set);
    // fprintf(stderr, "cpu=%d\n", cpu);
#endif
  }

  /**
   * @brief Execute a parallel for loop
   * @param ntasks Number of tasks to execute
   * @param f Function to execute for each task (takes task index as parameter)
   * @param thread_backend Threading backend to use (FTK_THREAD_*)
   * @param nthreads Number of threads (default: hardware_concurrency)
   * @param affinity Whether to set CPU affinity (default: false)
   */
  static void parallel_for(int ntasks, std::function<void(int)> f,
      int thread_backend = FTK_THREAD_PTHREAD,
      int nthreads = std::thread::hardware_concurrency(),
      bool affinity = false)
  {
    if (thread_backend == FTK_THREAD_PTHREAD) {
      nthreads = std::min(ntasks, nthreads);

      std::vector<std::thread> workers;
      for (auto i = 1; i < nthreads; i ++) {
        workers.push_back(std::thread([=]() {
          if (affinity) set_affinity(i);
          for (auto j = i; j < ntasks; j += nthreads)
            f(j);
        }));
      }

      if (affinity) set_affinity(0); // FIXME: for MPI runs, it seems that all work is congested on cpu0
      for (auto j = 0; j < ntasks; j += nthreads) // the main thread
        f(j);

      std::for_each(workers.begin(), workers.end(), [](std::thread &t) {t.join();});
    } else if (thread_backend == FTK_THREAD_OPENMP) {
#if FTK_HAVE_OPENMP
      fprintf(stderr, "parallelization w/ openmp...\n");
#pragma omp parallel for
      for (size_t j = 0; j < ntasks; j ++)
        f(j);
#else
      ftk::fatal(FTK_ERR_NOT_BUILT_WITH_OPENMP);
#endif
    } else if (thread_backend == FTK_THREAD_TBB) {
#if FTK_HAVE_TBB
      fprintf(stderr, "executing parallel_for with tbb...\n");
      tbb::parallel_for(tbb::blocked_range<size_t>(0, ntasks),
          [=](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++ i) 
              f(i);
          });
#else
      ftk::fatal(FTK_ERR_NOT_BUILT_WITH_TBB);
#endif
    } else 
      ftk::fatal(FTK_ERR_ACCELERATOR_UNSUPPORTED);
  }

  template <typename T, typename Container=std::set<T>>
  static void parallel_for(const Container& set, std::function<void(const T&)> f, 
      int xl, int nthreads, bool affinity) {
    std::vector<T> vector(set.size());
    std::copy(set.begin(), set.end(), vector.begin());

    parallel_for(set.size(), [&](int i) { f(vector[i]); }, 
        xl, nthreads, affinity);
  }

  template <typename Container> // =std::map<K, V>>
  static void parallel_for_container(Container& map, std::function<void(typename Container::iterator)> f,
      int xl = FTK_THREAD_PTHREAD, 
      int nthreads = std::thread::hardware_concurrency(), 
      bool affinity = true)
  {
    std::vector<typename Container::iterator> its;
    for (typename Container::iterator it = map.begin(); it != map.end(); it ++)
      its.push_back(it);

    parallel_for(its.size(), [&](int i) { f(its[i]); }, 
        xl, nthreads, affinity);
  }
 
protected:
  diy::mpi::communicator comm;
  int root_proc = 0;
};

}

#endif
