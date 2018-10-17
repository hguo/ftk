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
