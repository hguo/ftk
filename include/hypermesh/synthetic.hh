#ifndef _HYPERMESH_SYNTHETIC_DATA_H
#define _HYPERMESH_SYNTHETIC_DATA_H

namespace hypermesh {

// the synthetic woven function
template <typename T> 
T woven_function_2Dt(T x, T y, T t) 
{
  return cos(x*cos(t)-y*sin(t))*sin(x*sin(t)+y*cos(t));
}

// generate 2D woven data
template <typename T>
hypermesh::ndarray<T> synthetic_woven_2D(int DW, int DH, T t = T(1e-4), T scaling_factor = T(15))
{
  hypermesh::ndarray<T> scalar;
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

// generate 2D time-varying woven data
template <typename T>
hypermesh::ndarray<T> synthetic_woven_2Dt(int DW, int DH, int DT, T scaling_factor = T(15))
{
  hypermesh::ndarray<T> scalar;
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

}

#endif
