#ifndef _FTK_OBJECT_HH
#define _FTK_OBJECT_HH

#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <set>
#include <functional>
#include <algorithm>
#include <ftk/ftk_config.hh>
#include <ftk/external/diy/mpi.hpp>

namespace ftk {

struct object {
  static void fatal(const std::string& str) {
    std::cerr << "FATAL: " << str << std::endl;
    exit(1);
  }

  static void warn(const std::string& str) {
    std::cerr << "WARN: " << str << std::endl;
  }

  static void parallel_for(int ntasks, int nthreads, std::function<void(int)> f) {
    nthreads = std::min(ntasks, nthreads);

    std::vector<std::thread> workers;
    for (auto i = 1; i < nthreads; i ++) {
      workers.push_back(std::thread([=]() {
        for (auto j = i; j < ntasks; j += nthreads)
          f(j);
      }));
    }

    for (auto j = 0; j < ntasks; j += nthreads) // the main thread
      f(j);

    std::for_each(workers.begin(), workers.end(), [](std::thread &t) {t.join();});
  }

  template <typename T>
  static void parallel_for(const std::set<T>& set, int nthreads, std::function<void(const T&)> f) {
    std::vector<T> vector(set.size());
    std::copy(set.begin(), set.end(), vector.begin());

    parallel_for(set.size(), nthreads, [&](int i) { f(vector[i]); } );
  }
};

}

#endif
