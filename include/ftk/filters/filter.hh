#ifndef _FTK_FILTER_HH
#define _FTK_FILTER_HH

#include <ftk/ftk_config.hh>
#include <ftk/external/diy/mpi.hpp>
#include <ftk/external/cxxopts.hpp>
#include <mutex>
#include <thread>
#include <cassert>

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
    nthreads = default_nthreads();
  }

  filter(int argc, char **argv) {
    parse_arguments(argc, argv);
  }

  virtual void update() = 0;

  void use_accelerator(int i) {xl = i;}

  void parse_arguments(int argc, char **argv) {
    std::string str_xl;

    cxxopts::Options options(argv[0]);
    options.add_options()
      ("nthreads", "number of threads", 
        cxxopts::value<int>(nthreads)->default_value(std::to_string(default_nthreads())))
      ("x,accelerator", "use accelerator: none|cuda|kokkos|openmp|sycl|tbb", 
        cxxopts::value<std::string>(str_xl)->default_value("none"));
    auto results = options.parse(argc, argv);

    if (str_xl == "none") {
    } else if (str_xl == "cuda") {
      xl = FTK_XL_CUDA;
    } else {
      fprintf(stderr, "[FTK] fatal: unknow/unsupported accelerator %s\n", str_xl.c_str());
      assert(false);
    }
  }

  int default_nthreads() const {
    if (comm.size() > 1) return 1; // use 1 thread per proc for mpi runs
    else return std::thread::hardware_concurrency(); 
  }

protected:
  diy::mpi::communicator comm;

  int xl = FTK_XL_NONE;
  int nthreads = 1;
  std::mutex mutex;
};

}

#endif
