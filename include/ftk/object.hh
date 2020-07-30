#ifndef _FTK_OBJECT_HH
#define _FTK_OBJECT_HH

#include <iostream>
#include <string>
#include <thread>
#include <vector>

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
    for (size_t i = 1; i < nthreads; i ++) {
      workers.push_back(std::thread([&]() {
        for (size_t j = i; j < ntasks; j += nthreads)
          f(j);
      }));
    }

    for (size_t j = 0; j < ntasks; j += nthreads) // the main thread
      f(j);

    std::for_each(workers.begin(), workers.end(), [](std::thread &t) {t.join();});
  }
};

}

#endif
