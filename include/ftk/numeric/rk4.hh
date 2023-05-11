#ifndef _FTK_RK4_HH
#define _FTK_RK4_HH

#include <ftk/config.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/quaternion.hh>
#include <functional>

namespace ftk {

template <typename T=double>
void angular_stepping(
    const T *x, 
    const T* v, // assuming v is tangential to the sphere
    const T h, T *xn)
{
  const T R = vector_2norm<3, T>(x); // radius
  
  T axis[3]; // rotation axis
  cross_product(x, v, axis);

  const T vm = vector_2norm<3>(v); // velocity magnitude
  const T w = vm / R; // angular velocity
  
  // fprintf(stderr, "x=%f, %f, %f, axis=%f, %f, %f, w=%f\n", 
  //     x[0], x[1], x[2], axis[0], axis[1], axis[2], w);
  // std::cerr << w << ", " << x[0] << ", " << x[1] << ", " << x[2] << std::endl;

  const T dw = h * w; // angular step
  axis_rotate_vector(axis, dw, x, xn);

  xn[3] = xn[3] + h;
}

template <typename T=double>
void angular_and_vertical_stepping( // composition of angular and vertical velocities
    const T *x, 
    const T* v, // assuming v[0], v[1], v[2] is tangential to the sphere
                // 5 channels: x[3] time, x[4] vertical vel
    const T h, T *x1)
{
  const T R = vector_2norm<3, T>(x); // radius
  
  T axis[3]; // rotation axis
  cross_product(x, v, axis);

  const T vm = vector_2norm<3>(v); // velocity magnitude
  const T w = vm / R; // angular velocity
  const T dw = h * w; // angular step

  axis_rotate_vector(axis, dw, x, x1);

  const T R1 = R + v[4] * h; // new radius
  for (int i = 0; i < 3; i ++)
    x1[i] = x1[i] * R1 / R;

  x1[3] = x[3] + h * v[3]; // new time
}


template <typename T=double>
void spherical_stepping(const T *x, const T *v, const T h, T *xn)
{
  const T R = vector_2norm<3, T>(x); // radius
  for (int k = 0; k < 3; k ++)
    xn[k] = x[k] + h * v[k];
  
  const T Rn = vector_2norm<3, T>(xn); // new radius
  for (int k = 0; k < 3; k ++)
    xn[k] = xn[k] / Rn * R;
}

template <typename T=double> 
bool spherical_rk1(T *x, std::function<bool(const T*, T*)> f, T h, T *v0 = nullptr)
{
  T v[4]; // velocity
  if (!f(x, v)) return false;

  if (v0) 
    for (int k = 0; k < 3; k ++)
      v0[k] = v[k];

  // spherical_stepping(x, v, h, x);
  angular_stepping(x, v, h, x);

  return true;
}

template <typename T=double> 
bool spherical_rk1_with_vertical_velocity(T *x, std::function<bool(const T*, T*)> f, T h, T *v0 = nullptr)
{
  T v[10]; // velocity with the 4th component vertical
  if (!f(x, v)) return false;

  if (v0) 
    for (int k = 0; k < 5; k ++)
      v0[k] = v[k];

  // spherical_stepping(x, v, h, x);
  angular_stepping(x, v, h, x);
  
  x[4] += v[4] * h; // the vertical component
  return true;
}

template <typename T=double>
bool spherical_rk4_with_vertical_velocity(T *x, std::function<bool(const T*, T*)> f, T h, T *v0 = nullptr)
{
  constexpr int nch = 5;
  T v[10]; // velocity with the 4th component vertical
  T k1[nch], k2[nch], k3[nch], k4[nch];

  const T x0[4] = {x[0], x[1], x[2], x[3]};

  if (!f(x, v)) return false;
  if (v0) 
    for (int k = 0; k < 10; k ++)
      v0[k] = v[k];

  // fprintf(stderr, "x0=%f, %f, %f, %f, v=%f, %f, %f, %f, %f\n", 
  //     x[0], x[1], x[2], x[3], 
  //     v[0], v[1], v[2], v[3], v[4]);

  // k1 = f(x);
  for (int k = 0; k < nch; k ++)
    k1[k] = v[k];

  // k2 = f(x + 0.5*h*k1);
  T x2[4];
  angular_and_vertical_stepping<T>(x0, k1, h/2, x2);
  if (!f(x2, k2)) return false;

  // k3 = f(x + 0.5*h*k2);
  T x3[4];
  angular_and_vertical_stepping<T>(x0, k2, h/2, x3);
  if (!f(x3, k3)) return false;

  // k4 = f(x + k3);
  T x4[4];
  angular_and_vertical_stepping<T>(x0, k3, h, x4);
  if (!f(x4, k4)) return false;

  T w[nch];
  for (int i = 0; i < nch; i ++)
    w[i] = (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6.0;

  angular_and_vertical_stepping<T>(x0, w, h, x);
 
#if 0
  fprintf(stderr, "x0=%f, %f, %f, %f, v=%f, %f, %f, %f, %f, x=%f, %f, %f, %f\n", 
      x0[0], x0[1], x0[2], x0[3], 
      v[0], v[1], v[2], v[3], v[4], 
      x[0], x[1], x[2], x[3]);
#endif
  return true;
}

template <typename T=double>
bool rk4(int nd, T *pt, std::function<bool(const T*, T*)> f, T h, T *v0 = nullptr) // optional output of v at v0
{
  T p0[nd]; 
  memcpy(p0, pt, sizeof(T)*nd); 
  
  T v[nd]; 

  // 1st rk step
  if (!f(pt, v)) return false; 
  T k1[nd]; 
  for (int i=0; i<nd; i++) k1[i] = h*v[i];

  if (v0)
    memcpy(v0, v, sizeof(T)*nd);
  
  // 2nd rk step
  for (int i=0; i<nd; i++) pt[i] = p0[i] + 0.5*k1[i]; 
  if (!f(pt, v)) return false; 
  T k2[nd]; 
  for (int i=0; i<nd; i++) k2[i] = h*v[i]; 

  // 3rd rk step
  for (int i=0; i<nd; i++) pt[i] = p0[i] + 0.5*k2[i]; 
  if (!f(pt, v)) return false; 
  T k3[nd]; 
  for (int i=0; i<nd; i++) k3[i] = h*v[i]; 

  // 4th rk step
  for (int i=0; i<nd; i++) pt[i] = p0[i] + k3[i]; 
  if (!f(pt, v)) return false; 
  for (int i=0; i<nd; i++) 
    pt[i] = p0[i] + (k1[i] + 2.0*(k2[i]+k3[i]) + h*v[i])/6.0; 

  return true; 
}

}

#endif
