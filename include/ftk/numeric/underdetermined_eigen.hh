#ifndef _FTK_UNDERDETERMINED_EIGEN_HH
#define _FTK_UNDERDETERMINED_EIGEN_HH

#include "polynomial.hh"
#include <algorithm>
#include <mps/mps.h>
#include <gmp.h>

namespace ftk {

void solve_polynomials(double * x, int n){
	mps_context * s = mps_context_new();
	mps_context_select_algorithm(s, MPS_ALGORITHM_SECULAR_GA);
	mps_monomial_poly * poly = mps_monomial_poly_new(s, n);
	for(int i=0; i<=n; i++){
		mps_monomial_poly_set_coefficient_d(s, poly, i, x[i], 0);
	}
	mps_context_set_input_poly (s, MPS_POLYNOMIAL (poly));
	mps_mpsolve (s);
	mps_output(s);
	if(poly) mps_polynomial_free(s, poly);
	mps_context_free(s);
}

// solve Ax + a = lambda * (Bx + b)
// return x[i] = P[i][lambda] / Q[lambda]
template <typename T>
bool underdetermined_eigen_3x3(const T A[3][3], const T a[3], const T B[3][3], const T b[3], T P[3][4], T Q[4]){
	// polynomials p[], p[i] is the coeff for x^i
	// (A - lambda*B)x = - (a - lambda*b) 
	// A_lambda * x = a_lambda
	T A_lambda[3][3][2];
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			A_lambda[i][j][0] = A[i][j];
			A_lambda[i][j][1] = -B[i][j];
		}
	}
	T a_lambda[3][2];
	for(int i=0; i<3; i++){
		a_lambda[i][0] = -a[i];
		a_lambda[i][1] = b[i];
	}
	// det_A = Q
	for(int i=0; i<=3; i++)
		Q[i] = 0;

	// intermediate result
	T tmp[4];
	T tmp2[4];
	// 012 + 
	polynomial_multiplication(A_lambda[0][0], 1, A_lambda[1][1], 1, tmp);
	polynomial_multiplication(tmp, 2, A_lambda[2][2], 1, tmp2);
	polynomial_addition_in_place(Q, 3, tmp2, 3);
	// 021 -
	polynomial_multiplication(A_lambda[0][0], 1, A_lambda[1][2], 1, tmp);
	polynomial_multiplication(tmp, 2, A_lambda[2][1], 1, tmp2);
	polynomial_subtraction_in_place(Q, 3, tmp2, 3);
	// 120 +
	polynomial_multiplication(A_lambda[0][1], 1, A_lambda[1][2], 1, tmp);
	polynomial_multiplication(tmp, 2, A_lambda[2][0], 1, tmp2);
	polynomial_addition_in_place(Q, 3, tmp2, 3);
	// 102 -
	polynomial_multiplication(A_lambda[0][1], 1, A_lambda[1][0], 1, tmp);
	polynomial_multiplication(tmp, 2, A_lambda[2][2], 1, tmp2);
	polynomial_subtraction_in_place(Q, 3, tmp2, 3);
	// 201 +
	polynomial_multiplication(A_lambda[0][2], 1, A_lambda[1][0], 1, tmp);
	polynomial_multiplication(tmp, 2, A_lambda[2][1], 1, tmp2);
	polynomial_addition_in_place(Q, 3, tmp2, 3);
	// 210 -
	polynomial_multiplication(A_lambda[0][2], 1, A_lambda[1][1], 1, tmp);
	polynomial_multiplication(tmp, 2, A_lambda[2][0], 1, tmp2);
	polynomial_subtraction_in_place(Q, 3, tmp2, 3);

	T adjugate_A[3][3][3] = {0};
	// 00 +
	polynomial_multiplication(A_lambda[1][1], 1, A_lambda[2][2], 1, tmp);
	polynomial_multiplication(A_lambda[2][1], 1, A_lambda[1][2], 1, tmp2);
	polynomial_subtraction(tmp, 2, tmp2, 2, adjugate_A[0][0]);
	// 01 = C10 = -1 * det
	polynomial_multiplication(A_lambda[2][1], 1, A_lambda[0][2], 1, tmp);
	polynomial_multiplication(A_lambda[0][1], 1, A_lambda[2][2], 1, tmp2);
	polynomial_subtraction(tmp, 2, tmp2, 2, adjugate_A[0][1]);
	// 02 +
	polynomial_multiplication(A_lambda[0][1], 1, A_lambda[1][2], 1, tmp);
	polynomial_multiplication(A_lambda[1][1], 1, A_lambda[0][2], 1, tmp2);
	polynomial_subtraction(tmp, 2, tmp2, 2, adjugate_A[0][2]);
	// 10 -
	polynomial_multiplication(A_lambda[2][0], 1, A_lambda[1][2], 1, tmp);
	polynomial_multiplication(A_lambda[1][0], 1, A_lambda[2][2], 1, tmp2);
	polynomial_subtraction(tmp, 2, tmp2, 2, adjugate_A[1][0]);
	// 11 +
	polynomial_multiplication(A_lambda[0][0], 1, A_lambda[2][2], 1, tmp);
	polynomial_multiplication(A_lambda[2][0], 1, A_lambda[0][2], 1, tmp2);
	polynomial_subtraction(tmp, 2, tmp2, 2, adjugate_A[1][1]);
	// 12 -
	polynomial_multiplication(A_lambda[1][0], 1, A_lambda[0][2], 1, tmp);
	polynomial_multiplication(A_lambda[0][0], 1, A_lambda[1][2], 1, tmp2);
	polynomial_subtraction(tmp, 2, tmp2, 2, adjugate_A[1][2]);
	// 20 +
	polynomial_multiplication(A_lambda[1][0], 1, A_lambda[2][1], 1, tmp);
	polynomial_multiplication(A_lambda[2][0], 1, A_lambda[1][1], 1, tmp2);
	polynomial_subtraction(tmp, 2, tmp2, 2, adjugate_A[2][0]);
	// 21 -
	polynomial_multiplication(A_lambda[2][0], 1, A_lambda[0][1], 1, tmp);
	polynomial_multiplication(A_lambda[0][0], 1, A_lambda[2][1], 1, tmp2);
	polynomial_subtraction(tmp, 2, tmp2, 2, adjugate_A[2][1]);
	// 22 +
	polynomial_multiplication(A_lambda[0][0], 1, A_lambda[1][1], 1, tmp);
	polynomial_multiplication(A_lambda[1][0], 1, A_lambda[0][1], 1, tmp2);
	polynomial_subtraction(tmp, 2, tmp2, 2, adjugate_A[2][2]);

	for(int i=0; i<3; i++){
		for(int j=0; j<=3; j++){
			P[i][j] = 0;
		}
	}
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			polynomial_multiplication(adjugate_A[i][j], 2, a_lambda[j], 1, tmp);
			polynomial_addition_in_place(P[i], 3, tmp, 3);
		}
	}
	return true;
}

// solve underdetermined_eigen_3x3 with constrain x[0]*x[1] = x[2]
// return result[i] = x[0], x[1], x[2], and lambda[i]
template <typename T>
int solve_underdetermined_eigen_3x3_with_constrain(const T A[3][3], const T a[3], const T B[3][3], const T b[3], T lambda[6], T result[6][3],
	const T epsilon = std::numeric_limits<T>::epsilon()){
	T P[3][4] = {0};
	T Q[4] = {0};
	underdetermined_eigen_3x3(A, a, B, b, P, Q);
	std::cout << std::endl;
	std::cout << std::endl;
	for(int i=0; i<3; i++){
		for(int j=0; j<=3; j++){
			std::cout << P[i][j] << " ";
		}
		std::cout << std::endl;	
	}

	std::cout << std::endl;
	for(int j=0; j<=3; j++){
		std::cout << Q[j] << " ";
	}
	std::cout << std::endl;	

	T poly[7];
	polynomial_multiplication(P[0], 3, P[1], 3, poly);
	T tmp[7];
	polynomial_multiplication(P[2], 3, Q, 3, tmp);
	polynomial_subtraction_in_place(poly, 6, tmp, 6);
	double root_real[6], root_im[6];
	solve_polynomials(poly, 6, root_real, root_im);
	// validate real root
	int count = 0;
	for(int i=0; i<6; i++){
		std::cout << root_real[i] << " " << root_im[i] << std::endl;
		if(fabs(root_im[i]) < epsilon){
			// real root
			T Qx = polynomial_evaluate(Q, 3, root_real[i]);
			lambda[count] = root_real[i];
			result[count][0] = polynomial_evaluate(P[0], 3, root_real[i]) / Qx;
			result[count][1] = polynomial_evaluate(P[1], 3, root_real[i]) / Qx;
			result[count][2] = polynomial_evaluate(P[2], 3, root_real[i]) / Qx;
			count ++;
		}
	}
	return count;
}

}
#endif