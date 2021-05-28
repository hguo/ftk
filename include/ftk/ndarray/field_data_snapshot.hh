#ifndef _FTK_FIELD_DATA_SNAPSHOT_HH
#define _FTK_FIELD_DATA_SNAPSHOT_HH

#include <ftk/config.hh>
#include <ftk/ndarray/ndarray_base.hh>

namespace ftk {
 
// field data snapshot for various different trackers. 
// every field is optional
struct field_data_snapshot_t {
  std::shared_ptr<ndarray_base> 
    scalar,  // scalar field
    vector,  // vector field
    wector,  // another vector field
    vorticity, // vorticity field
    jacobian, // jacobian field
    scalar1, // additional scalar field
    uv,  // two scalar field or complex field
    rho, // magnetude of complex field
    phi, // phase of complex field
    re,  // real part of complex field
    im;  // imag part of complex field
};

}

#endif
