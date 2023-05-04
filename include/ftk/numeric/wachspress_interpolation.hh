#ifndef _FTK_WACHSPRESS_INTERPOLATION_HH
#define _FTK_WACHSPRESS_INTERPOLATION_HH

#include <ftk/config.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>

namespace ftk {

template <typename T=double>
__device__ __host__
inline T triangle_area3(const T x0[3], const T x1[3], const T x2[3])
{
  const T u[3] = {x1[0] - x0[0], x1[1] - x0[1], x1[2] - x0[2]}, 
          v[3] = {x2[0] - x0[0], x2[1] - x0[1], x2[2] - x0[2]};

  T w[3];
  cross_product(u, v, w);
  
  return 0.5 * vector_2norm<3, T>(w);
}

template <typename T=double>
__device__ __host__ 
inline bool wachspress_weights( // return true if x is inside the cell
    const int nverts,
    const T X[][3], 
    const T x[3], // cartisian coordinates
    T omega[]) // wachspress coordinates
{
  T A[nverts]; // areas of triangles that connect x
  for (int i = 0; i < nverts; i ++)
  A[i] = triangle_area3(X[i], x, X[(i+1)%nverts]);
   
  // fprintf(stderr, "A=%f, %f, %f, %f, %f, %f\n",
  //     A[0], A[1], A[2], A[3], A[4], A[5]);

  for (int i = 0; i < nverts; i ++) {
    const T B = triangle_area3(X[i], X[(i-1+nverts)%nverts], X[(i+1)%nverts]);

    T prodA(1); // summation of areas on the opposite side
    for (int j = 0; j < nverts; j ++) {
      // if (i == j || i == (j-1+nverts)%nverts)
      if (j == i || j == (i-1+nverts)%nverts)
        continue;
      else
        prodA *= A[j];
    }

    omega[i] = B * prodA;
    
    // fprintf(stderr, "--B%d=%f, prodA=%f\n", i, B, prodA);
  }

  T sum_omega(0);
  for (int i = 0; i < nverts; i ++)
    sum_omega += omega[i];
  
  // fprintf(stderr, "omega=%f, %f, %f, %f, %f, %f, sum=%f\n", 
  //     omega[0], omega[1], omega[2], omega[3], omega[4], omega[5], sum_omega);
  
  for (int i = 0; i < nverts; i ++)
    omega[i] = omega[i] / sum_omega;
  
  // fprintf(stderr, "omega=%f, %f, %f, %f, %f, %f\n", 
  //     omega[0], omega[1], omega[2], omega[3], omega[4], omega[5]);
    
  return true;
}

#if 0
template <typename T=double, int nvars=3>
__device__ __host__ 
inline bool wachspress_interpolation3(
    const int nverts,
    const T X[][3], const T V[][nvars], // inputs
    const T x[3], // coordinates
    T f[nvars]) // outputs
{
  return false;
}
#endif

}

#endif
