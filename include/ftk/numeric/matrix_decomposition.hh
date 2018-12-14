#ifndef _FTK_MATRIX_DECOMPOSITION_HH
#define _FTK_MATRIX_DECOMPOSITION_HH

namespace ftk {

template <class T>
inline void matrix_symmetric_decomposition_3x3(const T M[3][3], T S[3][3])
{
  S[0][0] = M[0][0]; 
  S[0][1] = T(0.5) * (M[0][1] + M[1][0]);
  S[0][2] = T(0.5) * (M[0][2] + M[2][0]);
  S[1][0] = T(0.5) * (M[1][0] + M[0][1]);
  S[1][1] = M[1][1];
  S[1][2] = T(0.5) * (M[1][2] + M[2][1]);
  S[2][0] = T(0.5) * (M[2][0] + M[0][2]);
  S[2][1] = T(0.5) * (M[2][1] + M[1][2]);
  S[2][2] = M[2][2];
}

template <class T>
inline void matrix_antisymmetric_decomposition_3x3(const T M[3][3], T A[3][3])
{
  A[0][0] = T(0);
  A[0][1] = T(0.5) * (M[0][1] - M[1][0]);
  A[0][2] = T(0.5) * (M[0][2] - M[2][0]);
  A[1][0] = T(0.5) * (M[1][0] - M[0][1]);
  A[1][1] = T(0);
  A[1][2] = T(0.5) * (M[1][2] - M[2][1]);
  A[2][0] = T(0.5) * (M[2][0] - M[0][2]);
  A[2][1] = T(0.5) * (M[2][1] - M[1][2]);
  A[2][2] = T(0);
}

}

#endif
