#ifndef _TDGL_META_HH
#define _TDGL_META_HH

#include <ftk/config.hh>

namespace ftk {

struct tdgl_metadata_t {
  int ndims; 
  int dims[3];
  bool pbc[3];
  float zaniso;
  float lengths[3], origins[3], cell_lengths[3];
  float time;
  float B[3];
  float Jxext, Kex, Kex_dot, V;
  float fluctuation_amp;
  int dtype;
}; 

}

#endif
