#ifndef _FTK_RBF_HH
#define _FTK_RBF_HH

#include <ftk/config.hh>
#include <ftk/numeric/vector_dot_product.hh>
#include <ftk/numeric/vector_norm.hh>

namespace ftk {

template <typename T=double>
__device__ __host__
inline void rbf3d_plane_vec_const_dir(
    const int n,
    const T X[][3], // src points
    const T N[][3], // (unit) normal vectors
    const T x[3],   // dst point
    const T alpha,  // radius
    const T t[2][3],   // planar basis vector of the tangent plane
    T coeffs[][3])
{
  const int m = n + 2;

  T planarSourcePoints[n][2], planarUnitVectors[n][2];
  for (int i = 0; i < n; i ++) {
    planarSourcePoints[i][0] = vector_dot_product3(X[i], t[0]);
    planarSourcePoints[i][1] = vector_dot_product3(X[i], t[1]);
    planarUnitVectors[i][0] = vector_dot_product3(N[i], t[0]);
    planarUnitVectors[i][1] = vector_dot_product3(N[i], t[1]);
  }

  const T planarDestinationPoint[2] = {
    vector_dot_product3(x, t[0]), 
    vector_dot_product3(x, t[1])
  };

  T matrix[m*m];
  const T alpha2 = alpha * alpha;
  for (int i = 0; i < n; i ++) {
    for (int j = i; j < n; j ++) {
      const T r2 = vector_dist_2norm2<2, T>(planarSourcePoints[i], planarSourcePoints[j]) / alpha2;
      const T rbf = 1.0 / (r2 + 1);
      const T dot = vector_dot_product2(planarUnitVectors[i], planarUnitVectors[j]);

      matrix[i*m+j] = rbf * dot;
      matrix[j*m+i] = matrix[i*m+j];
    }
  }
  for (int i = 0; i < n; i ++) {
    for (int j = 0; j < 2; j ++) {
      matrix[i*m + n+j] = planarUnitVectors[i][j];
      matrix[(n+j)*m + i] = matrix[i*m + n+j];
    }
  }

  T rhs[2][m];
  for (int i = 0; i < n; i ++) {
    const T r2 = vector_dist_2norm2<2, T>(planarDestinationPoint, planarSourcePoints[i]) / alpha2;
    for (int j = 0; j < 2; j ++)
      rhs[j][i] = 1.0 / (r2 + 1) * planarUnitVectors[i][j];
  }
  for (int i = 0; i < 2; i ++)
    rhs[i][n+i] = 1.0;

  T coeffs0[2][m];
  legsc(m, matrix, rhs[0], coeffs0[0]);
  legsc(m, matrix, rhs[1], coeffs0[1]);

  for (int i = 0; i < 3; i ++) 
    for (int j = 0; j < n; j ++)
      coeffs[j][i] = t[0][i] * coeffs0[0][j] + t[1][i] * coeffs0[1][j];
}

}

#endif
