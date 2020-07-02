#ifndef _FTK_INVERSE_TRILINEAR_INTERPOLATION_SOLVER_HH
#define _FTK_INVERSE_TRILINEAR_INTERPOLATION_SOLVER_HH

#include <ftk/numeric/underdetermined_eigen_solver.hh>

namespace ftk {

// find the zero point in [0, 1]x[0, 1]x[0, 1] quad, using generalized eigenvalue problem
template <typename T>
int inverse_trilinear_interpolation8_3(const T V[8][3], T result[6][3]){
	const T *V000 = V[0], *V001 = V[1], *V010 = V[2], *V011 = V[3],
		  *V100 = V[4], *V101 = V[5], *V110 = V[6], *V111 = V[7];

	// Axyz + Bxy + Cxz + Dyz + Ex + Fy + Gz + H = 0
	T A[3], B[3], C[3], D[3], E[3], F[3], G[3], H[3];
	for(int i=0; i<3; i++){
		A[i] = -V000[i] + V001[i] + V010[i] - V011[i] + V100[i] - V101[i] - V110[i] + V111[i];
		B[i] = V000[i] - V010[i] - V100[i] + V110[i];
		C[i] = V000[i] - V001[i] - V100[i] + V101[i];
		D[i] = V000[i] - V001[i] - V010[i] + V011[i];
		E[i] = -V000[i] + V100[i];
		F[i] = -V000[i] + V010[i];
		G[i] = -V000[i] + V001[i];
		H[i] = V000[i];
	}
	// transform to A'x+a' = lambda*(B'x + b')
	// where x = (x, y, xy), lambda = z
	T Ap[3][3], ap[3], Bp[3][3], bp[3];
	for(int i=0; i<3; i++) Ap[i][0] = E[i];
	for(int i=0; i<3; i++) Ap[i][1] = F[i];
	for(int i=0; i<3; i++) Ap[i][2] = B[i];
	for(int i=0; i<3; i++) ap[i] = H[i];
	for(int i=0; i<3; i++) Bp[i][0] = - C[i];
	for(int i=0; i<3; i++) Bp[i][1] = - D[i];
	for(int i=0; i<3; i++) Bp[i][2] = - A[i];
	for(int i=0; i<3; i++) bp[i] = - G[i];
	T lambda[6];
	int root_count = solve_underdetermined_eigen_3x3_with_constrain(Ap, ap, Bp, bp, lambda, result);
	return root_count;
}

}
#endif
