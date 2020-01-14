#ifndef _HYPERMESH_SYNTHETIC_DATA_H
#define _HYPERMESH_SYNTHETIC_DATA_H

#include <ftk/ndarray.hh>
#include <ftk/hypermesh/lattice.hh>

namespace ftk {

// the synthetic woven function
template <typename T> 
T woven_function_2Dt(T x, T y, T t) 
{
  return cos(x*cos(t)-y*sin(t))*sin(x*sin(t)+y*cos(t));
}

// generate 2D woven data
template <typename T>
ndarray<T> synthetic_woven_2D(int DW, int DH, T t = T(1e-4), T scaling_factor = T(15))
{
  ndarray<T> scalar;
  scalar.reshape(DW, DH);

  // const T scaling_factor = 15; // the factor that controls the shape of the synthesize data
  for (int j = 0; j < DH; j ++) {
    for (int i = 0; i < DW; i ++) {
      const T x = ((T(i) / (DW-1)) - 0.5) * scaling_factor,
              y = ((T(j) / (DH-1)) - 0.5) * scaling_factor;
      scalar(i, j) = woven_function_2Dt(x, y, t);
    }
  }

  return scalar;
}

// generate 2D woven data in a core
template <typename T>
ndarray<T> synthetic_woven_2D_part(const lattice& ext, const lattice& core, T t = T(1e-4), T scaling_factor = T(15))
{
  ndarray<T> scalar;
  scalar.reshape(core.sizes());

  const int DW = ext.size(0), DH = ext.size(1);
  for (int j = 0; j < core.size(1) ; j ++) {
    for (int i = 0; i < core.size(0); i ++) {
      const T x = ((T(i + core.start(0)) / (DW-1)) - 0.5) * scaling_factor,
              y = ((T(j + core.start(1)) / (DH-1)) - 0.5) * scaling_factor;
      scalar(i, j) = woven_function_2Dt(x, y, t);
    }
  }

  return scalar;
}

// generate 2D time-varying woven data
template <typename T>
ndarray<T> synthetic_woven_2Dt(int DW, int DH, int DT, T scaling_factor = T(15))
{
  ndarray<T> scalar;
  scalar.reshape(DW, DH, DT);

  // const T scaling_factor = 15; // the factor that controls the shape of the synthesize data
  for (int k = 0; k < DT; k ++) {
    for (int j = 0; j < DH; j ++) {
      for (int i = 0; i < DW; i ++) {
        const T x = ((T(i) / (DW-1)) - 0.5) * scaling_factor,
                y = ((T(j) / (DH-1)) - 0.5) * scaling_factor, 
                t = (T(k) / (DT-1)) + 1e-4;
        scalar(i, j, k) = woven_function_2Dt(x, y, t);
      }
    }
  }

  return scalar;
}

// double gyre 2D flow
template <typename T>
ndarray<T> synthetic_double_gyre(int DW, int DH)
{
  // TODO
}

// ABC flow
template <typename T>
ndarray<T> synthetic_abc_flow(int DW, int DH, int DD, 
    T A=std::sqrt(T(3)), T B=std::sqrt(T(2)), T C=T(1))
{

}

// 2D merger
template <typename T>
T merger_function_2Dt(T x, T y, T t)
{
  auto f = [](T cx, T cy, T x, T y) {return exp(-((x-cx)*(x-cx) + (y-cy)*(y-cy)));};

  // add rotation
  T xp = x * cos(t) - y * sin(t), 
    yp = x * sin(t) + y * cos(t);
  x = xp;
  y = yp;

  T cx0 = sin(t - M_PI_2), // + 1e-4, 
    cx1 = sin(t + M_PI_2), // + 1e-4,
    cy0 = 1e-4,
    cy1 = 1e-4;
        
  return std::max(f(cx0, cy0, x, y), f(cx1, cy1, x, y));
}

template <typename T>
ndarray<T> synthetic_merger_2D(int DW, int DH, T t)
{
  ndarray<T> scalar;
  scalar.reshape(DW, DH);
  
  for (int j = 0; j < DH; j ++) {
    for (int i = 0; i < DW; i ++) {
      // the domain is [-2, 2]x[-2, 2]
      const T x = ((T(i) / (DW-1)) - 0.5) * 4,
              y = ((T(j) / (DH-1)) - 0.5) * 4;
      scalar(i, j) = merger_function_2Dt(x, y, t);
    }
  }

  return scalar;
}

}

#endif
