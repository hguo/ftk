#ifndef _FTK_MPAS_PARTICLE_HH
#define _FTK_MPAS_PARTICLE_HH

#include <ftk/config.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/quaternion.hh>

namespace ftk {

template <typename T=double>
struct mpas_particle_t {
  T r[3]; // unit vector pointing from the center towards the surface
  T R; // distance from the center
  T t; // time
  T v[3]; // velocity
  T vv;   // vertical velocity
  T scalar[FTK_CP_MAX_NUM_VARS]; // scalar attributes
  int hint_c = 0, hint_l = 0;

  T get_z() const { return R - R0; }
  void set_z(T z) { R = R0 + z; }

  void get_x(T x[3]) const {
    for (int i = 0; i < 3; i ++)
      x[i] = r[i] * R;
  }
  void set_x(const T x[3]) {
    R = vector_2norm<3>(x);
    for (int i = 0; i < 3; i ++)
      r[i] = x[i] / R;
  }

  void axis_rotate(const T axis[3], const T rad) { // axis rotation
    T rn[3]; // new orientation
    axis_rotate_vector(axis, rad, r, rn);
    for (int i = 0; i < 3; i ++)
      r[i] = rn[i];
  }

  void step(const T dt)
  {
    T axis[3]; // rotation axis
    cross_product(r, v, axis);

    const T vm = vector_2norm<3>(v); // velocity magnitude
    const T w = vm / R; // angular velocity
    const T theta = dt * w;
    
    axis_rotate(axis, theta);
    R += vv * dt;
    t += dt;
  }

  static constexpr T R0 = 6371229.0;
};

}

#endif
