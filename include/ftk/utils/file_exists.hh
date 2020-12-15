#ifndef _FTK_FILE_EXISTS
#define _FTK_FILE_EXISTS

#include <iostream>
#include <string>

namespace ftk {

bool file_exists(const std::string& filename) {
  std::ifstream f(filename);
  return f.good();
}

bool file_not_exists(const std::string& filename) { 
  return !file_exists(filename); 
}

}

#endif
