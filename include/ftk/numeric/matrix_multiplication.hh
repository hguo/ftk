#ifndef _FTK_MATRIX_MULTIPLICATION_HH
#define _FTK_MATRIX_MULTIPLICATION_HH

namespace ftk {

template <class T, int M, int K, int N>
void matrix_matrix_multiplication(const T A[M][K], const T B[K][N], T C[M][N])
{
  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++) {
      T dot = T(0);
      for (int k=0; k<K; k++)
        dot += A[i][k] * B[k][j];
      C[i][j] = dot;
    }
  }
}

template <class T>
void matrix3x3_scalar_multiplication(const T A[9], T b, T C[9])
{
  for (int i=0; i<9; i++) 
    C[i] = A[i]*b;
}

template <class T>
void matrix2x2_vector2_multiplication(const T A[2][2], const T b[2], T x[2])
{
  x[0] = A[0][0] * b[0] + A[0][1] * b[1];
  x[1] = A[1][0] * b[0] + A[1][1] * b[1];
}

template <class T>
void matrix2x3_vector3_multiplication(const T A[2][3], const T b[3], T x[2])
{
  x[0] = A[0][0] * b[0] + A[0][1] * b[1] + A[0][2] * b[2];
  x[1] = A[1][0] * b[0] + A[1][1] * b[1] + A[1][2] * b[2];
}

template <class T>
__device__ __host__
void matrix3x3_vector3_multiplication(const T A[9], const T b[3], T c[3]) 
{
  c[0] = A[0] * b[0] + A[1] * b[1] + A[2] * b[2];
  c[1] = A[3] * b[0] + A[4] * b[1] + A[5] * b[2];
  c[2] = A[6] * b[0] + A[7] * b[1] + A[8] * b[2];
}

template <class T>
__device__ __host__
void matrix3x3_vector3_multiplication(const T A[3][3], const T b[3], T c[3])
{
  c[0] = A[0][0] * b[0] + A[0][1] * b[1] + A[0][2] * b[2];
  c[1] = A[1][0] * b[0] + A[1][1] * b[1] + A[1][2] * b[2];
  c[2] = A[2][0] * b[0] + A[2][1] * b[1] + A[2][2] * b[2];
}

template <class T>
void matrix2x2_matrix2x2_multiplication(const T A[2][2], const T B[2][2], T C[2][2])
{
  C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];
  C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1];
  C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];
  C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1];
}

template <class T>
void matrix2x2_matrix2x3_multiplication(const T A[2][2], const T B[2][3], T C [2][3])
{
  C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];
  C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1];
  C[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2];
  C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];
  C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1];
  C[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2];
}

template <class T>
void matrix3x3_matrix3x3_multiplication(const T A[9], const T B[9], T C[9])
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      T dot = T(0);
      for (int k=0; k<3; k++)
        dot += A[i*3+k] * B[k*3+j];
      C[i*3+j] = dot;
    }
  }
}

template <class T>
void matrix3x3_matrix3x3_multiplication(const T A[3][3], const T B[3][3], T C[3][3])
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      T dot = T(0);
      for (int k=0; k<3; k++)
        dot += A[i][k] * B[k][j];
      C[i][j] = dot;
    }
  }
}

template <class T>
inline void matrix3x3_matrix3x4_multiplication(const T A[3][3], const T B[3][4], T C[3][4])
{
  matrix_matrix_multiplication<T, 3, 3, 4>(A, B, C);
}

template <class T>
inline void matrix3x3_square(const T M[3][3], T M2[3][3])
{
  matrix_matrix_multiplication_3x3_3x3(M, M, M2);
}

template <class T>
inline void matrix2x3_matrix3x2_multiplication(const T A[2][3], const T B[3][2], T C[2][2])
{
  C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
  C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1];
  C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0];
  C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1];
}

template <class T>
inline void matrix3x4_matrix4x3_multiplication(const T A[3][4], const T B[4][3], T C[3][3])
{
  matrix_matrix_multiplication<T, 3, 4, 3>(A, B, C);
}

template <class T>
inline void matrix3x4_matrix4x2_multiplication(const T A[3][4], const T B[4][2], T C[2][2])
{
  matrix_matrix_multiplication<T, 3, 4, 2>(A, B, C);
}

}

#endif
