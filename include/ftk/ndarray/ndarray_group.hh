#ifndef _FTK_NDARRAY_GROUP_HH
#define _FTK_NDARRAY_GROUP_HH

#include <ftk/ndarray.hh>
#include <unordered_map>

namespace ftk {

struct ndarray_group : public std::unordered_map<std::string, ndarray_base*> {
  ndarray_group() {}

  // template <typename ... Args> ndarray_group(Args&&... args);
};

}

#endif
