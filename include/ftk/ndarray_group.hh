#ifndef _FTK_NDARRAY_GROUP_HH
#define _FTK_NDARRAY_GROUP_HH

#include <ftk/ndarray.hh>

namespace ftk {

struct ndarray_group : public std::map<std::string, ndarray_base> {
  ndarray_group() {}

  // template <typename ... Args> ndarray_group(Args&&... args);
};

}

#endif
