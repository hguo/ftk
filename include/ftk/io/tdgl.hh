#ifndef _TDGL_IO_HH
#define _TDGL_IO_HH

#include <ftk/ftk_config.hh>

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
  tdgl_reader(const std::string& filename, 
      bool derive_rho_phi_ = true, 
      bool derive_re_im_ = true,
      bool derive_J_ = true) :
    derive_rho_phi(derive_rho_phi_), 
    derive_re_im(derive_re_im_),
    derive_J(derive_J_) {}

  bool read();

public:
  bool derive_rho_phi, derive_re_im, derive_J;

  tdgl_metadata_t meta;
  ndarray<float> rho_phi, re_im, J;
};

}

#endif
