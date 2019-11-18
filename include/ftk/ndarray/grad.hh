#ifndef _HYPERMESH_GRAD_HH
#define _HYPERMESH_GRAD_HH

#include <ftk/ndarray.hh>

namespace ftk {

// derive 2D gradients for 2D scalar field
template <typename T>
ndarray<T> gradient2D(const ndarray<T>& scalar)
{
  const int DW = scalar.dim(0), DH = scalar.dim(1);
  ndarray<T> grad;
  grad.reshape(2, DW, DH); 
  
  for (int j = 1; j < DH-1; j ++) {
    for (int i = 1; i < DW-1; i ++) {
      auto dfdx = grad(0, i, j) = 0.5 * (scalar(i+1, j) - scalar(i-1, j)) * (DW-1);
      auto dfdy = grad(1, i, j) = 0.5 * (scalar(i, j+1) - scalar(i, j-1)) * (DH-1);
      // fprintf(stderr, "s=%f, grad=%f, %f\n", scalar(i, j), dfdx, dfdy);
    }
  }
  return grad;
}

// derive 2D gradients for 2D time varying scalar field
template <typename T>
ndarray<T> gradient2Dt(const ndarray<T>& scalar)
{
  const int DW = scalar.dim(0), DH = scalar.dim(1), DT = scalar.dim(2);
  ndarray<T> grad;
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
ndarray<T> jacobian2D(const ndarray<T>& vec)
{
  const int DW = vec.dim(1), DH = vec.dim(2);
  ndarray<T> grad;
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
ndarray<T> jacobian2Dt(const ndarray<T>& vec)
{
  const int DW = vec.dim(1), DH = vec.dim(2), DT = vec.dim(3);
  ndarray<T> grad;
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

// derive gradients for 3D scalar field
template <typename T>
ndarray<T> gradient3D(const ndarray<T>& scalar)
{
  const int DW = scalar.dim(0), DH = scalar.dim(1), DD = scalar.dim(2);
  ndarray<T> grad;
  grad.reshape(3, DW, DH, DD);

  for (int k = 1; k < DD-1; k ++) {
    for (int j = 1; j < DH-1; j ++) {
      for (int i = 1; i < DW-1; i ++) {
        grad(0, i, j, k) = 0.5 * (scalar(i+1, j, k) - scalar(i-1, j, k));
        grad(1, i, j, k) = 0.5 * (scalar(i, j+1, k) - scalar(i, j-1, k));
        grad(2, i, j, k) = 0.5 * (scalar(i, j, k+1) - scalar(i, j, k-1));
      }
    }
  }

  return grad;
}

template <typename T>
ndarray<T> gradient3Dt(const ndarray<T>& scalar)
{
  const int DW = scalar.dim(0), DH = scalar.dim(1), DD = scalar.dim(2), DT = scalar.dim(3);
  ndarray<T> grad;
  grad.reshape(3, DW, DH, DD, DT);

  for (int t = 0; t < DT; t ++) {
    for (int k = 1; k < DD-1; k ++) {
      for (int j = 1; j < DH-1; j ++) {
        for (int i = 1; i < DW-1; i ++) {
          grad(0, i, j, k, t) = 0.5 * (scalar(i+1, j, k, t) - scalar(i-1, j, k, t));
          grad(1, i, j, k, t) = 0.5 * (scalar(i, j+1, k, t) - scalar(i, j-1, k, t));
          grad(2, i, j, k, t) = 0.5 * (scalar(i, j, k+1, t) - scalar(i, j, k-1, t));
        }
      }
    }
  }

  return grad;
}

// derivate gradients (jacobians) for 3D vector field
template <typename T>
ndarray<T> jacobian3D(const ndarray<T>& V)
{
  const int DW = V.dim(1), DH = V.dim(2), DD = V.dim(3);
  ndarray<T> J;
  J.reshape({3, 3, DW, DH, DD});

  for (int k = 0; k < DD-2; k ++) {
    for (int j = 2; j < DH-2; j ++) {
      for (int i = 2; i < DW-2; i ++) {
        const float H00 = J(0, 0, i, j, k) = // ddf/dx2
          0.5 * (V(0, i+1, j, k) - V(0, i-1, j, k));
        const float H01 = J(0, 1, i, j, k) = // ddf/dxdy
          0.5 * (V(0, i, j+1, k) - V(0, i, j-1, k));
        const float H02 = J(0, 2, i, j, k) = // ddf/dxdz
          0.5 * (V(0, i, j, k+1) - V(0, i, j, k-1));

        const float H10 = J(1, 0, i, j, k) = // ddf/dydx
          0.5 * (V(1, i+1, j, k) - V(1, i-1, j, k));
        const float H11 = J(1, 1, i, j, k) = // ddf/dy2
          0.5 * (V(1, i, j+1, k) - V(1, i, j-1, k));
        const float H12 = J(1, 2, i, j, k) = // ddf/dydz
          0.5 * (V(1, i, j, k+1) - V(1, i, j, k-1));

        const float H20 = J(2, 0, i, j, k) = // ddf/dydx
          0.5 * (V(2, i+1, j, k) - V(2, i-1, j, k));
        const float H21 = J(2, 1, i, j, k) = // ddf/dy2
          0.5 * (V(2, i, j+1, k) - V(2, i, j-1, k));
        const float H22 = J(2, 2, i, j, k) = // ddf/dydz
          0.5 * (V(2, i, j, k+1) - V(2, i, j, k-1));
      }
    }
  }

  return J;
}

// derive gradients (jacobians) for 3D time varying vector field
template <typename T>
ndarray<T> jacobian3Dt(const ndarray<T>& V)
{
  const int DW = V.dim(1), DH = V.dim(2), DD = V.dim(3), DT = V.dim(4);
  ndarray<T> J;
  J.reshape(3, 3, DW, DH, DD, DT);

  for (int t = 0; t < DT; t ++) {
    for (int k = 0; k < DD-2; k ++) {
      for (int j = 2; j < DH-2; j ++) {
        for (int i = 2; i < DW-2; i ++) {
          const float H00 = J(0, 0, i, j, k, t) = // ddf/dx2
            0.5 * (V(0, i+1, j, k, t) - V(0, i-1, j, k, t));
          const float H01 = J(0, 1, i, j, k) = // ddf/dxdy
            0.5 * (V(0, i, j+1, k, t) - V(0, i, j-1, k, t));
          const float H02 = J(0, 2, i, j, k) = // ddf/dxdz
            0.5 * (V(0, i, j, k+1, t) - V(0, i, j, k-1, t));

          const float H10 = J(1, 0, i, j, k, t) = // ddf/dydx
            0.5 * (V(1, i+1, j, k, t) - V(1, i-1, j, k, t));
          const float H11 = J(1, 1, i, j, k, t) = // ddf/dy2
            0.5 * (V(1, i, j+1, k, t) - V(1, i, j-1, k, t));
          const float H12 = J(1, 2, i, j, k, t) = // ddf/dydz
            0.5 * (V(1, i, j, k+1, t) - V(1, i, j, k-1, t));

          const float H20 = J(2, 0, i, j, k, t) = // ddf/dydx
            0.5 * (V(2, i+1, j, k, t) - V(2, i-1, j, k, t));
          const float H21 = J(2, 1, i, j, k, t) = // ddf/dy2
            0.5 * (V(2, i, j+1, k, t) - V(2, i, j-1, k, t));
          const float H22 = J(2, 2, i, j, k, t) = // ddf/dydz
            0.5 * (V(2, i, j, k+1, t) - V(2, i, j, k-1, t));
        }
      }
    }
  }

  return J;
}

}

#endif
