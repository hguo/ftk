#ifndef _FTK_INVERSE_BILINEAR_INTERPOLATION_SOLVER_HH
#define _FTK_INVERSE_BILINEAR_INTERPOLATION_SOLVER_HH

#include <ftk/ftk_config.hh>

namespace ftk {

// find the zero point in [0, 1]x[0, 1] quad, using generalized eigenvalue problem
template <typename T>
bool inverse_bilinear_interpolation(const T u[4], const T v[4], T pos[2])
{
  // find_zero_unit_quad_bilinear(const T re[4], const T im[4], T pos[2], T epsilon=0)
  T f00 = u[0], f10 = u[1], f01 = u[3], f11 = u[2], // counter-clockwise
    g00 = v[0], g10 = v[1], g01 = v[3], g11 = v[2];
  T A0 = f00 - f10 - f01 + f11,  // Axy + Bx + Cy + D = 0
    B0 = f10 - f00, 
    C0 = f01 - f00, 
    D0 = f00,
    A1 = g00 - g10 - g01 + g11, 
    B1 = g10 - g00, 
    C1 = g01 - g00, 
    D1 = g00; 
  T M0[4] = {-B0, -D0, -B1, -D1}, // stored in row major
    M1[4] = {A0, C0, A1, C1}; // (yM1 - M0)v = 0, v = {x, 1}^T

  T detM1 = A0*C1 - A1*C0; // TODO: check if detM1==0
  T invM1[4] = {C1/detM1, -C0/detM1, -A1/detM1, A0/detM1};  

  // Q = invM1*M0
  T Q[4] = {
    invM1[0]*M0[0] + invM1[1]*M0[2], 
    invM1[0]*M0[1] + invM1[1]*M0[3], 
    invM1[2]*M0[0] + invM1[3]*M0[2], 
    invM1[2]*M0[1] + invM1[3]*M0[3]
  };

  // compute y=eig(Q)
  T trace = Q[0] + Q[3];
  T det = Q[0]*Q[3] - Q[1]*Q[2];
  T lambda[2] = {
    static_cast<T>(trace/2 + sqrt(trace*trace/4 - det)), 
    static_cast<T>(trace/2 - sqrt(trace*trace/4 - det))
  }; 

  T x[2] = {
    (lambda[0]-Q[3])/Q[2], 
    (lambda[1]-Q[3])/Q[2]
  }; 
  T y[2] = {
    lambda[0], 
    lambda[1]
  };

  T xx, yy;
  bool found = false; 
  for (int i=0; i<2; i++) // check the two roots 
    if (x[i]>=0 && x[i]<=1 && y[i]>=0 && y[i]<=1) {
      pos[0] = x[i]; 
      pos[1] = y[i];
      found = true; 
      break; 
    }
  
  return found; 
}

template <typename T>
bool inserve_barycentric3_linear_interpolation3(T v[3][3], T w[3][3], T lambda[4])
{
  return false;
}

} // namespace ftk

#endif
