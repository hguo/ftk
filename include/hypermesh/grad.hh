#ifndef _HYPERMESH_GRAD_HH
#define _HYPERMESH_GRAD_HH

#include "ndarray.hh"

namespace hypermesh {

// derive 2D gradients for 2D scalar field
template <typename T>
hypermesh::ndarray<T> gradient2D(const hypermesh::ndarray<T>& scalar)
{
  const int DW = scalar.dim(0), DH = scalar.dim(1);
  hypermesh::ndarray<T> grad;
  grad.reshape(2, DW, DH); 
  
  for (int j = 1; j < DH-1; j ++) {
    for (int i = 1; i < DW-1; i ++) {
      grad(0, i, j) = 0.5 * (scalar(i+1, j) - scalar(i-1, j)) * (DW-1);
      grad(1, i, j) = 0.5 * (scalar(i, j+1) - scalar(i, j-1)) * (DH-1);
    }
  }
  return grad;
}

// derive 2D gradients for 2D time varying scalar field
template <typename T>
hypermesh::ndarray<T> gradient2Dt(const hypermesh::ndarray<T>& scalar)
{
  const int DW = scalar.dim(0), DH = scalar.dim(1), DT = scalar.dim(2);
  hypermesh::ndarray<T> grad;
  grad.reshape(2, DW, DH, DT);
  
  for (int k = 0; k < DT; k ++) {
    for (int j = 1; j < DH-1; j ++) {
      for (int i = 1; i < DW-1; i ++) {
        grad(0, i, j, k) = 0.5 * (scalar(i+1, j, k) - scalar(i-1, j, k)) * (DW-1);
        grad(1, i, j, k) = 0.5 * (scalar(i, j+1, k) - scalar(i, j-1, k)) * (DH-1);
      }
    }
  }
  return grad;
}

// derive gradients for 2D vector field
template <typename T>
hypermesh::ndarray<T> jacobian2D(const hypermesh::ndarray<T>& vec)
{
  const int DW = vec.dim(1), DH = vec.dim(2);
  hypermesh::ndarray<T> grad;
  grad.reshape(2, 2, DW, DH);

  for (int j = 2; j < DH-2; j ++) {
    for (int i = 2; i < DW-2; i ++) {
      const T H00 = grad(0, 0, i, j) = // du/dx 
        0.5 * (vec(0, i+1, j) - vec(0, i-1, j)) * (DW-1); 
      const T H01 = grad(0, 1, i, j) = // du/dy
        0.5 * (vec(0, i, j+1) - vec(0, i, j-1)) * (DH-1);
      const T H10 = grad(1, 0, i, j) = // dv/dx
        0.5 * (vec(1, i+1, j) - vec(1, i-1, j)) * (DW-1);
      const T H11 = grad(1, 1, i, j) = // dv.dy
        0.5 * (vec(1, i, j+1) - vec(1, i, j-1)) * (DH-1);
    }
  }
  return grad;
}

// derive 2D gradients for 2D vector field
template <typename T>
hypermesh::ndarray<T> jacobian2Dt(const hypermesh::ndarray<T>& vec)
{
  const int DW = vec.dim(1), DH = vec.dim(2), DT = vec.dim(3);
  hypermesh::ndarray<T> grad;
  grad.reshape(2, 2, DW, DH, DT);

  for (int k = 0; k < DT; k ++) {
    for (int j = 2; j < DH-2; j ++) {
      for (int i = 2; i < DW-2; i ++) {
        const T H00 = grad(0, 0, i, j, k) = // ddf/dx2
          0.5 * (vec(0, i+1, j, k) - vec(0, i-1, j, k)) * (DW-1); 
        const T H01 = grad(0, 1, i, j, k) = // ddf/dxdy
          0.5 * (vec(0, i, j+1, k) - vec(0, i, j-1, k)) * (DH-1);
        const T H10 = grad(1, 0, i, j, k) = // ddf/dydx
          0.5 * (vec(1, i+1, j, k) - vec(1, i-1, j, k)) * (DW-1);
        const T H11 = grad(1, 1, i, j, k) = // ddf/dy2
          0.5 * (vec(1, i, j+1, k) - vec(1, i, j-1, k)) * (DH-1);
      }
    }
  }
  return grad;
}

}

#endif
