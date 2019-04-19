#ifndef _HYPERMESH_CONV_HH
#define _HYPERMESH_CONV_HH

#include <cmath>
#include <string>
#include "ndarray.hh"

namespace hypermesh {

template <typename T>
ndarray<T> gaussian_kernel2D(T sigma, size_t ksizex, size_t ksizey)
{
  ndarray<T> kernel({ksizex, ksizey});

  // fill the kernel
  T centerx = static_cast<double>(ksizex - 1) * .5,
         centery = static_cast<double>(ksizey - 1) * .5;
  T r, s = T(2) * sigma * sigma;
  T sum(0);

  for (auto j = 0; j < ksizey; ++j) {
    for (auto i = 0; i < ksizex; ++i) {
      T x = static_cast<T>(i) - centerx,
        y = static_cast<T>(j) - centery;
      r = std::sqrt(x * x + y * y);
      // kernel[j * ksizex + i] = (std::exp(-(r * r) / s)) / (M_PI * s);
      kernel(i, j) = (std::exp(-(r * r) / s)) / (M_PI * s);
      sum += kernel[j * ksizex + i];
    }
  }

  // normalize the kernel
  T tmp = T(1) / sum;
  for (auto i = 0; i < ksizex * ksizey; ++i)
    kernel[i] *= tmp;

  return kernel;
}

template <typename T>
ndarray<T> conv2D(const ndarray<T> &data, const ndarray<T> &kernel, size_t padding = 0)
{
  const auto dimx = data.dim(0), 
             dimy = data.dim(1);
  const auto ksizex = kernel.dim(0), 
             ksizey = kernel.dim(1);
  // dimensions for the resulting data
  const auto dimx_r = dimx + padding * 2 - ksizex + 1,
             dimy_r = dimy + padding * 2 - ksizey + 1;

  // resulting data
  ndarray<T> res({dimx_r, dimy_r});
  // T* res = new T[dimx_r * dimy_r];

  // convolution
  // uncomment the following line to enable openmp
#pragma omp parallel for collapse(2)
  for (auto y = 0; y < dimy_r; ++y) {
    for (auto x = 0; x < dimx_r; ++x) {

      for (auto ky = 0; ky < ksizey; ++ky) {
        for (auto kx = 0; kx < ksizex; ++kx) {
          auto realy = y - padding + ky,
              realx = x - padding + kx;
          if (realx >= 0 && realx < dimx && realy >= 0 && realy < dimy) {
            res(x, y) += data(realx, realy) * kernel(kx, ky);
            // res[y * dimx_r + x] += data[realy * dimx + realx] *
            //                        kernel[ky * ksizex + kx];
          }
        }
      }

      res(x, y) /= ksizey * ksizex;
      // res[y * dimx_r + x] /= ksizey * ksizex;
    }
  }

  return res;
}

template <typename T>
ndarray<T> conv2D_gaussian(const ndarray<T> &data, 
    T sigma, size_t ksizex = 5, size_t ksizey = 5, size_t padding = 0)
{
  // make kernel
  const auto kernel = gaussian_kernel2D(sigma, ksizex, ksizey);
  // T* kernel = new T[ksizex * ksizey];
  // MakeGaussianKernel2D(sigma, ksizex, ksizey, kernel);

  // convolution
  auto res = conv2D(data, kernel, padding);
  // T* res = Conv2D(data, dimx, dimy, kernel, ksizex, ksizey, padding);

  // delete kernel
  // delete[] kernel;

  return res;
}

template <typename T>
T* Conv3D(T* data, int dimx, int dimy, int dimz,
          T* kernel, int ksizex, int ksizey, int ksizez,
          int padding=0) {
  // dimensions for the resulting data
  int dimx_r = dimx + padding * 2 - ksizex + 1,
      dimy_r = dimy + padding * 2 - ksizey + 1,
      dimz_r = dimz + padding * 2 - ksizez + 1;

  // resulting data
  T* res = new T[dimx_r * dimy_r * dimz_r];

  // convolution
  // uncomment the following line to enable openmp
  // #pragma omp parallel for collapse(3)
  for (int z = 0; z < dimz_r; ++z) {
    for (int y = 0; y < dimy_r; ++y) {
      for (int x = 0; x < dimx_r; ++x) {

        for (int kz = 0; kz < ksizez; ++kz) {
          for (int ky = 0; ky < ksizey; ++ky) {
            for (int kx = 0; kx < ksizex; ++kx) {
              int realz = z - padding + kz,
                  realy = y - padding + ky,
                  realx = x - padding + kx;
              if (realx >= 0 && realx < dimx &&
                  realy >= 0 && realy < dimy &&
                  realz >= 0 && realz < dimz) {
                res[z * dimy_r * dimx_r + y * dimx_r + x] +=
                    data[realz * dimy * dimx + realy * dimx + realx] *
                    kernel[kz * ksizey * ksizex + ky * ksizex + kx];
              }
            }
          }
        }

        res[z * dimy_r * dimx_r + y * dimx_r + x] /=
            ksizez * ksizey * ksizex;
      }
    }
  }

  return res;
}

}  // namespace hypermesh

#endif  // _HYPERMESH_CONV_HH
