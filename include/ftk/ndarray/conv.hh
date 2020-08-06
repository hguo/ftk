#ifndef _HYPERMESH_CONV_HH
#define _HYPERMESH_CONV_HH

#include <cmath>
#include <string>
#include <ftk/ndarray.hh>

namespace ftk {

// 2D convolutions
template <typename T>
ndarray<T> conv2D(const ndarray<T> &data, const ndarray<T> &kernel,
                  size_t padding = 0)
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
          }
        }
      }

      res(x, y) /= ksizey * ksizex;
    }
  }

  return res;
}

template <typename T>
std::vector<T> gaussian_kernel(T sigma, size_t size)
{
  std::vector<T> kernel(size);

  double center = static_cast<double>(size - 1) * 0.5;
  double r, s = 2. * sigma * sigma;
  double sum = 0.;

  for (auto i = 0; i < size; ++i) {
    double x = static_cast<double>(i) - center;
    r = x * x;
    kernel[i] = std::exp(-r / s);
    sum += kernel[i];
  }

  for (auto i = 0; i < size; i ++) 
    kernel[i] /= sum;

  // for (auto i = 0; i < size; i ++)
  //   fprintf(stderr, "kernel[%d]=%f\n", i, kernel[i]);

  return kernel;
}

template <typename T>
ndarray<T> gaussian_kernel2D(T sigma, size_t ksizex, size_t ksizey)
{
  ndarray<T> kernel({ksizex, ksizey});

  // fill the kernel
  double centerx = static_cast<double>(ksizex - 1) * .5,
         centery = static_cast<double>(ksizey - 1) * .5;
  double r, s = 2. * sigma * sigma;
  double sum = 0.;

  for (auto j = 0; j < ksizey; ++j) {
    for (auto i = 0; i < ksizex; ++i) {
      double x = static_cast<double>(i) - centerx,
             y = static_cast<double>(j) - centery;
      r = x * x + y * y;
      kernel(i, j) = std::exp(-r / s);
      sum += kernel(i, j);
    }
  }

  // normalize the kernel
  for (auto i = 0; i < ksizex * ksizey; ++i)
    kernel[i] /= sum;

  return kernel;
}

template <typename T>
ndarray<T> conv2D_gaussian(const ndarray<T> &data, T sigma,
                           size_t ksizex = 5, size_t ksizey = 5,
                           size_t padding = 0)
{
  // make kernel
  const auto kernel = gaussian_kernel2D(sigma, ksizex, ksizey);

  // convolution
  auto res = conv2D(data, kernel, padding);

  return res;
}

// 3D convolutions
template <typename T>
ndarray<T> conv3D(const ndarray<T> &data, const ndarray<T> &kernel,
                  size_t padding = 0)
{
  const auto dimx = data.dim(0),
             dimy = data.dim(1),
             dimz = data.dim(2);
  const auto ksizex = kernel.dim(0),
             ksizey = kernel.dim(1),
             ksizez = kernel.dim(2);
  // dimensions for the resulting data
  const auto dimx_r = dimx + padding * 2 - ksizex + 1,
             dimy_r = dimy + padding * 2 - ksizey + 1,
             dimz_r = dimz + padding * 2 - ksizez + 1;

  // resulting data
  ndarray<T> res({dimx_r, dimy_r, dimz_r});

  // convolution
  // uncomment the following line to enable openmp
#pragma omp parallel for collapse(3)
  for (auto z = 0; z < dimz_r; ++z) {
    for (auto y = 0; y < dimy_r; ++y) {
      for (auto x = 0; x < dimx_r; ++x) {

        for (auto kz = 0; kz < ksizez; ++kz) {
          for (auto ky = 0; ky < ksizey; ++ky) {
            for (auto kx = 0; kx < ksizex; ++kx) {
              auto realz = z - padding + kz,
                   realy = y - padding + ky,
                   realx = x - padding + kx;
              if (realx >= 0 && realx < dimx &&
                  realy >= 0 && realy < dimy &&
                  realz >= 0 && realz < dimz) {
                res(x, y, z) += data(realx, realy, realz) * kernel(kx, ky, kz);
              }
            }
          }
        }

        res(x, y, z) /= ksizez * ksizey * ksizex;
      }
    }
  }

  return res;
}

template <typename T>
ndarray<T> gaussian_kernel3D(
    T sigma, size_t ksizex, size_t ksizey, size_t ksizez)
{
  ndarray<T> kernel({ksizex, ksizey, ksizez});

  // fill the kernel
  double centerx = static_cast<double>(ksizex - 1) * .5,
         centery = static_cast<double>(ksizey - 1) * .5,
         centerz = static_cast<double>(ksizez - 1) * .5;
  double r, s = 2. * sigma * sigma;
  double sum = 0.;

  for (int j = 0; j < ksizey; ++j) {
    for (int i = 0; i < ksizex; ++i) {
      for (int k = 0; k < ksizez; ++k) {
        double x = static_cast<double>(i) - centerx,
               y = static_cast<double>(j) - centery,
               z = static_cast<double>(k) - centerz;
        r = x * x + y * y + z * z;
        kernel(i, j, k) = std::exp(-r/ s);
        sum += kernel(i, j, k);
      }
    }
  }

  // normalize the kernel
  for (auto i = 0; i < ksizex * ksizey* ksizez; ++i)
    kernel[i] /= sum;

  return kernel;
}

template <typename T>
ndarray<T> conv3D_gaussian(
    const ndarray<T> &data, T sigma,
    size_t ksizex = 5, size_t ksizey = 5, size_t ksizez = 5,
    size_t padding = 0)
{
  // make kernel
  const auto kernel = gaussian_kernel3D(sigma, ksizex, ksizey, ksizez);

  // convolution
  auto res = conv3D(data, kernel, padding);

  return res;
}

template <typename T>
ndarray<T> conv_gaussian(
    const ndarray<T> &data, T sigma,
    size_t ksize=5, size_t padding=0)
{
  if (data.nd() == 2) return conv2D_gaussian<T>(data, sigma, ksize, ksize, padding);
  else if (data.nd() == 3) return conv3D_gaussian<T>(data, sigma, ksize, ksize, ksize, padding);
  else return ndarray<T>();
}

}  // namespace ftk

#endif  // _HYPERMESH_CONV_HH
