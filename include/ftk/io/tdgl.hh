#ifndef _TDGL_IO_HH
#define _TDGL_IO_HH

#include <ftk/config.hh>
#include <ftk/ndarray.hh>
#include <ftk/io/tdgl_metadata.hh>
#include <iostream>

namespace ftk {

struct tdgl_reader {
  tdgl_reader(const std::string& filename_, 
      bool read_field_data_ = true,
      bool derive_J_ = false) : 
    filename(filename_),
    read_field_data(read_field_data_),
    derive_J(derive_J_) {}

  bool read();

  const tdgl_metadata_t& get_meta() const {return meta;}

public:
  bool read_field_data, derive_J;

  const std::string& filename;
  tdgl_metadata_t meta;
  ndarray<float> rho, phi, re, im, J;
};

}

#endif
