#ifndef _TDGL_IO_HH
#define _TDGL_IO_HH

#include <ftk/ftk_config.hh>
#include <ftk/ndarray.hh>
#include <iostream>

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

struct tdgl_reader {
  tdgl_reader(const std::string& filename_, 
      bool read_field_data_ = true,
      bool derive_J_ = false) : 
    filename(filename_),
    read_field_data(read_field_data_),
    derive_J(derive_J_) {}

  bool read();

public:
  bool read_field_data, derive_J;

  const std::string& filename;
  tdgl_metadata_t meta;
  ndarray<float> rho, phi, re, im, J;
};

}

#endif
