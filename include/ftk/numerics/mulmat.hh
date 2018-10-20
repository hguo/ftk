#ifndef _FTK_MULMAT
#define _FTK_MULMAT

namespace ftk {

template <class ValueType>
void mulmat3s(const ValueType A[9], ValueType b, ValueType C[9])
{
  for (int i=0; i<9; i++) 
    C[i] = A[i]*b;
}

template <class ValueType>
void mulmat3v(const ValueType A[9], const ValueType b[3], ValueType c[3]) 
{
  c[0] = A[0] * b[0] + A[1] * b[1] + A[2] * b[2];
  c[1] = A[3] * b[0] + A[4] * b[1] + A[5] * b[2];
  c[2] = A[6] * b[0] + A[7] * b[1] + A[8] * b[2];
}

template <class ValueType>
void mulmat3(const ValueType A[9], const ValueType B[9], ValueType C[9])
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      ValueType dot = ValueType(0);
      for (int k=0; k<3; k++)
        dot += A[i*3+k] * B[k*3+j];
      C[i*3+j] = dot;
    }
  }
}

}

#endif
