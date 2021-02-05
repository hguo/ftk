#ifndef _FTK_NDARRAY_GROUP_HH
#define _FTK_NDARRAY_GROUP_HH

#include <ftk/ndarray.hh>

namespace ftk {

struct ndarray_group : std::map<std::string, ndarray_wrapper> {
};

}

#endif
