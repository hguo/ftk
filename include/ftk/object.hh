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

enum { // thread backend
  FTK_THREAD_NONE = 0,
  FTK_THREAD_PTHREAD = 0,
  FTK_THREAD_OPENMP = 1,
  FTK_THREAD_TBB = 6
};

enum { 
  FTK_XL_NONE = 0,
  FTK_XL_SYCL = 2,
  FTK_XL_CUDA = 4,
  FTK_XL_KOKKOS_CUDA = 5
};

struct object {
  object() {}
  object(diy::mpi::communicator c) {comm = c;}

  void set_communicator(const diy::mpi::communicator comm_) {comm = comm_;}
  void set_root_proc(int p) {root_proc = p;}
  int get_root_proc() const {return root_proc;}
  bool is_root_proc() const {return root_proc == comm.rank();}

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

  static void parallel_for(int ntasks, std::function<void(int)> f, 
      int thread_backend = FTK_THREAD_PTHREAD, 
      int nthreads = std::thread::hardware_concurrency(), 
      bool affinity = true)
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

      if (affinity) set_affinity(0);
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
      fatal(FTK_ERR_NOT_BUILT_WITH_OPENMP);
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
      fatal(FTK_ERR_NOT_BUILT_WITH_TBB);
#endif
    } else 
      fatal(FTK_ERR_ACCELERATOR_UNSUPPORTED);
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
