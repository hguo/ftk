#ifndef _HYPERMESH_CONV_HH
#define _HYPERMESH_CONV_HH

namespace hypermesh {

template <typename T>
T* Conv2D(T* data, int dimx, int dimy,
          T* kernel, int ksizex, int ksizey,
          int padding=0) {
  // dimensions for the resulting data
  int dimx_r = dimx + padding * 2 - ksizex + 1,
      dimy_r = dimy + padding * 2 - ksizey + 1;

  // resulting data
  T* res = new T[dimx_r * dimy_r];

  // convolution
  // uncomment the following line to enable openmp
  // #pragma omp parallel for collapse(2)
  for (int y = 0; y < dimy_r; ++y) {
    for (int x = 0; x < dimx_r; ++x) {

      for (int ky = 0; ky < ksizey; ++ky) {
        for (int kx = 0; kx < ksizex; ++kx) {
          int realy = y - padding + ky,
              realx = x - padding + kx;
          if (realx >= 0 && realx < dimx && realy >= 0 && realy < dimy) {
            res[y * dimx_r + x] += data[realy * dimx + realx] *
                                   kernel[ky * ksizex + kx];
          }
        }
      }

      res[y * dimx_r + x] /= ksizey * ksizex;
    }
  }

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
