#ifndef _FTK_MATRIX_MULTIPLICATION_HH
#define _FTK_MATRIX_MULTIPLICATION_HH

namespace ftk {

template <class T>
void matrix_scalar_multiplication_3x3(const T A[9], T b, T C[9])
{
  for (int i=0; i<9; i++) 
    C[i] = A[i]*b;
}

template <class T>
void matrix_vector_multiplication_3x3(const T A[9], const T b[3], T c[3]) 
{
  c[0] = A[0] * b[0] + A[1] * b[1] + A[2] * b[2];
  c[1] = A[3] * b[0] + A[4] * b[1] + A[5] * b[2];
  c[2] = A[6] * b[0] + A[7] * b[1] + A[8] * b[2];
}

template <class T>
void matrix_vector_multiplication_3x3(const T A[3][3], const T b[3], T c[3])
{
  c[0] = A[0][0] * b[0] + A[0][1] * b[1] + A[0][2] * b[2];
  c[1] = A[1][0] * b[0] + A[1][1] * b[1] + A[1][2] * b[2];
  c[2] = A[2][0] * b[0] + A[2][1] * b[1] + A[2][2] * b[2];
}

template <class T>
void matrix_matrix_multiplication_2x2_2x2(const T A[2][2], const T B[2][2], T C[2][2])
{
  C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];
  C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1];
  C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];
  C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1];
}

template <class T>
void matrix_matrix_multiplication_3x3_3x3(const T A[9], const T B[9], T C[9])
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
void matrix_matrix_multiplication_3x3_3x3(const T A[3][3], const T B[3][3], T C[3][3])
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
void matrix_square_3x3(const T M[3][3], T M2[3][3])
{
  matrix_matrix_multiplication_3x3_3x3(M, M, M2);
}

}

#endif
