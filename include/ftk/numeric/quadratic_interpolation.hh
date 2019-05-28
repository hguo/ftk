#ifndef _FTK_QUADRATIC_INTERPOLATION_HH
#define _FTK_QUADRATIC_INTERPOLATION_HH

namespace ftk {

template <typename T>
void quadratic_interpolation_coefficient_simplex2(const T f[6], T Q[3][3])
{
  Q[0][0] = f[0];
  Q[0][1] = 2 * f[3] - (f[0] + f[1]) / 2;
  Q[0][2] = 2 * f[5] - (f[0] + f[2]) / 2;
  Q[1][0] = Q[0][1];
  Q[1][1] = f[1];
  Q[1][2] = 2 * f[4] - (f[1] + f[2]) / 2;
  Q[2][0] = Q[0][2];
  Q[2][1] = Q[1][2];
  Q[2][2] = f[2];
}

template <typename T>
void quadratic_interpolation_gradient_simplex2(const T Q[3][3], const T mu[3], T grad[3])
{
  grad[0] = 2 * (Q[0][0]*mu[0] + Q[0][1]*mu[1] + Q[0][2]*mu[2]);
  grad[1] = 2 * (Q[1][1]*mu[1] + Q[0][1]*mu[0] + Q[1][2]*mu[2]);
  grad[2] = 2 * (Q[2][2]*mu[2] + Q[0][2]*mu[0] + Q[1][2]*mu[1]);
}

template <typename T>
void quadratic_interpolation_gradient_simplex2_vertices(const T Q[3][3], T grad[3][3])
{
  // gradient at 100
  grad[0][0] = 2 * Q[0][0];
  grad[0][1] = 2 * Q[0][1];
  grad[0][2] = 2 * Q[0][2];

  // gradient at 010
  grad[1][0] = 2 * Q[1][0];
  grad[1][1] = 2 * Q[1][1];
  grad[1][2] = 2 * Q[1][2];
  
  // gradient at 001
  grad[2][0] = 2 * Q[2][0];
  grad[2][1] = 2 * Q[2][1];
  grad[2][2] = 2 * Q[2][2];
}


// coefficient of quadratic interpolation
// f: scalar fields
// x, y: coordinates of x1 x2 x3
template <typename T>
bool quadratic_interpolation_coefficients(const T f[6], const T x[3][2], T Q[3][3]){
	/*
				x_3
			x_6		x_5
		x_1		x_4		x_2
	*/
	// det(x1x2x3)
	// TODO: check det = 0
	T det_123 = x[0][0]*x[1][1] + x[1][0]*x[2][1] + x[2][0]*x[0][1] - x[0][0]*x[2][1] - x[1][0]*x[0][1] - x[2][0]*x[1][1];
	T L1_x = (x[1][1] - x[2][1]) / det_123;
	T L1_y = (x[2][0] - x[1][0]) / det_123;
	T L1_c = (x[1][0]*x[2][1] - x[2][0]*x[1][1]) / det_123;
	T L2_x = (x[2][1] - x[0][1]) / det_123;
	T L2_y = (x[0][0] - x[2][0]) / det_123;
	T L2_c = (x[2][0]*x[0][1] - x[0][0]*x[2][1]) / det_123;
	T L3_x = (x[0][1] - x[1][1]) / det_123;
	T L3_y = (x[1][0] - x[0][0]) / det_123;
	T L3_c = (x[0][0]*x[1][1] - x[1][0]*x[0][1]) / det_123;
	// quardratic
	// f = L1*(2*L1 - 1)*f(x1) + L2*(2*L2 - 1)*f(x2) + L3*(2*L3 - 1)*f(x3) + 
	//		4*L1*L2*f(x4) + 4*L2*L3*f(x5) + 4*L1*L3*f(x6)
	// f = Ax^2 + By^2 + Cxy + Dx + Ey + F = x^T Q x
	Q[0][0] = 2 * L1_x * L1_x * f[0] + 2 * L2_x * L2_x * f[1] + 2 * L3_x * L3_x * f[2] +
				4 * L1_x * L2_x * f[3] + 4 * L2_x * L3_x * f[4] + 4 * L1_x * L3_x * f[5];
	Q[1][1] = 2 * L1_y * L1_y * f[0] + 2 * L2_y * L2_y * f[1] + 2 * L3_y * L3_y * f[2] +
				4 * L1_y * L2_y * f[3] + 4 * L2_y * L3_y * f[4] + 4 * L1_y * L3_y * f[5];
	Q[0][1] = Q[1][0] = (4 * L1_x * L1_y * f[0] + 4 * L2_x * L2_y * f[1] + 4 * L3_x * L3_y * f[2] +
				4 * (L1_x * L2_y + L1_y * L2_x) * f[3] + 4 * (L2_x * L3_y + L2_y * L3_x) * f[4] + 
				+ 4 * (L1_x * L3_y + L1_y * L3_x) * f[5]) / 2;
	Q[0][2] = Q[2][0] = ((4 * L1_x * L1_c - L1_x) * f[0] + (4 * L2_x * L2_c - L2_x) * f[1] + (4 * L3_x * L3_c - L3_x) * f[2] +
				4 * (L1_x * L2_c + L1_c * L2_x) * f[3] + 4 * (L2_x * L3_c + L2_c * L3_x) * f[4] +
				4 * (L1_x * L3_c + L1_c * L3_x) * f[5]) / 2;
	Q[1][2] = Q[2][1] = ((4 * L1_y * L1_c - L1_y) * f[0] + (4 * L2_y * L2_c - L2_y) * f[1] + (4 * L3_y * L3_c - L3_y) * f[2] +
				4 * (L1_y * L2_c + L1_c * L2_y) * f[3] + 4 * (L2_y * L3_c + L2_c * L3_y) * f[4] +
				4 * (L1_y * L3_c + L1_c * L3_y) * f[5]) / 2;
	Q[2][2] = (2 * L1_c * L1_c - L1_c) * f[0] + (2 * L2_c * L2_c - L2_c) * f[1] + (2 * L3_c * L3_c - L3_c) * f[2] +
				4 * L1_c * L2_c * f[3] + 4 * L2_c * L3_c * f[4] + 4 * L1_c * L3_c * f[5];
	return true;
}

template <typename T>
T quadratic_interpolation(const T Q[3][3], const T x[2])
{
  return Q[0][0] * x[0] * x[0] 
    + T(2) * Q[0][1] * x[0] * x[1] 
    + T(2) * Q[0][2] * x[0]
    + Q[1][1] * x[1] * x[1]
    + T(2) * Q[1][2] * x[1]
    + Q[2][2];
}

/**
 * @param grad  gradient vector on three vertices
 */
template <typename T>
bool quadratic_interpolation_gradient(const T Q[3][3], const T x[3][2], T grad[3][2]) {
	for(int i=0; i<3; i++){
		grad[i][0] = 2 * Q[0][0] * x[i][0] + (Q[0][1] + Q[1][0]) * x[i][1] + (Q[0][2] + Q[2][0]);
		grad[i][1] = (Q[0][1] + Q[1][0]) * x[i][0] + 2 * Q[1][1] * x[i][1] + (Q[1][2] + Q[2][1]);
	}
	return true;
}

template <typename T>
bool quadratic_interpolation_hessian(const T Q[3][3], T H[2][2]){
	H[0][0] = Q[0][0];
	H[1][1] = Q[1][1];
	H[0][1] = H[1][0] = Q[0][1] + Q[1][0];
	return true;
}

}

#endif
