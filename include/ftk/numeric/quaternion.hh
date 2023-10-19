#ifndef _FTK_QUATERNION_HH
#define _FTK_QUATERNION_HH

#include <ftk/config.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_dot_product.hh>
#include <ftk/numeric/vector_normalization.hh>
#include <cmath>

namespace ftk {

template <typename T>
__device__ __host__
void quatmul(const T x[], const T y[], T z[]) {
	const T ix[3] = {x[1], x[2], x[3]}, 
					iy[3] = {y[1], y[2], y[3]};
  T cross[3];
	cross_product(ix, iy, cross);
	const T dot = vector_dot_product3(ix, iy);
	
  z[0] = x[0] * y[0] - dot, 
	z[1] = x[0] * iy[0] + y[0] * ix[0] + cross[0],
	z[2] = x[0] * iy[1] + y[0] * ix[1] + cross[1],
	z[3] = x[0] * iy[2] + y[0] * ix[2] + cross[2];
}

template <typename T>
__device__ __host__
void axis_rotate_vector(
    const T axis_[3], 
    const T theta, 
    const T x[3], 
    T x1[3])
{
  T axis[3];
  vector_normalization2<3, T>(axis_, axis);

  const T c = std::cos(theta / 2), 
          s = std::sin(theta / 2);

  const T q[4] = {c, s * axis[0], s * axis[1], s * axis[2]};
  const T invq[4] = {q[0], -q[1], -q[2], -q[3]};
  const T p[4] = {0, x[0], x[1], x[2]};

  T q1[4], q2[4];
  quatmul(q, p, q1);
  quatmul(q1, invq, q2);

  x1[0] = q2[1]; 
  x1[1] = q2[2];
  x1[2] = q2[3];
}

}

#endif
