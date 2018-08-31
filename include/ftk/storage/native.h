#ifndef _FTK_DIR_STORAGE
#define _FTK_DIR_STORAGE

#include "ftk/storage/base.h"
#include <sys/stat.h>
#include <errno.h>

namespace ftk {

class storage_native: public storage {
public: 
  bool open(const std::string& dbname) {
    return false;
  }

  void close() {
  }

  void put(const std::string& key, const std::string& val) {
  }

  std::string get(const std::string& key) {
    std::string val;
    return val;
  }

private:
};

}

#endif
