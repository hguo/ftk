#ifndef _FTK_NDARRAY_GROUP_HH
#define _FTK_NDARRAY_GROUP_HH

#include <ftk/ndarray.hh>
#include <ftk/ndarray_wrapper.hh>

namespace ftk {

struct ndarray_group : public std::map<std::string, ndarray_wrapper> {
  ndarray_group() {}

  // template <typename ... Args> ndarray_group(Args&&... args);
};

}

#endif
