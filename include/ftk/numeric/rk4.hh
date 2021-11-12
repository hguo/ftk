#ifndef _FTK_RK4_HH
#define _FTK_RK4_HH

#include <functional>

namespace ftk {

template <int nd, typename T=double>
bool rk4(T *pt, std::function<bool(const T*, T*)> f, T h)
{
  T p0[nd]; 
  memcpy(p0, pt, sizeof(T)*nd); 
  
  T v[nd]; 

  // 1st rk step
  if (!f(pt, v)) return false; 
  T k1[nd]; 
  for (int i=0; i<nd; i++) k1[i] = h*v[i]; 
  
  // 2nd rk step
  for (int i=0; i<nd; i++) pt[i] = p0[i] + 0.5*k1[i]; 
  if (!f(pt, v)) return true; 
  T k2[nd]; 
  for (int i=0; i<nd; i++) k2[i] = h*v[i]; 

  // 3rd rk step
  for (int i=0; i<nd; i++) pt[i] = p0[i] + 0.5*k2[i]; 
  if (!f(pt, v)) return true; 
  T k3[nd]; 
  for (int i=0; i<nd; i++) k3[i] = h*v[i]; 

  // 4th rk step
  for (int i=0; i<nd; i++) pt[i] = p0[i] + k3[i]; 
  if (!f(pt, v)) return true; 
  for (int i=0; i<nd; i++) 
    pt[i] = p0[i] + (k1[i] + 2.0*(k2[i]+k3[i]) + h*v[i])/6.0; 

  return true; 
}

}

#endif
