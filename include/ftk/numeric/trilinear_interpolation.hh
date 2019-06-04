#ifndef _FTK_TRILINEAR_INTERPOLATION_H
#define _FTK_TRILINEAR_INTERPOLATION_H

namespace ftk {

template <typename T>
inline void trilinear_interpolation3(const T V[8][3], const T alpha, const T beta, const T gamma, T result[3]){
	for(int i=0; i<3; i++){
		// V000 = V[0], V001 = V[1], V010 = V[2], V011 = V[3]
		// V100 = V[4], V101 = V[5], V110 = V[6], V111 = V[7]
		T c00 = (1-gamma) * V[0][i] + gamma * V[1][i];
		T c01 = (1-gamma) * V[2][i] + gamma * V[3][i];
		T c10 = (1-gamma) * V[4][i] + gamma * V[5][i];
		T c11 = (1-gamma) * V[6][i] + gamma * V[7][i];
		T c0 = (1-beta) * c00 + beta * c01;
		T c1 = (1-beta) * c10 + beta * c11;
		result[i] = (1-alpha) * c0 + alpha * c1;
	}
}

}
#endif