#ifndef _FTK_TRANSPOSE_H
#define _FTK_TRANSPOSE_H

namespace ftk {

template <typename ValueType>
void transpose2(ValueType m[4])
{
  swap(m[1], m[2]);
}

template <typename ValueType>
void transpose3(ValueType m[9]) 
{
  swap(m[1], m[3]);
  swap(m[2], m[6]);
  swap(m[5], m[7]);
}

template <typename ValueType>
void transpose3(const ValueType a[9], ValueType b[9]) 
{
  b[0] = a[0];
  b[1] = a[3];
  b[2] = a[6];
  b[3] = a[1];
  b[4] = a[4];
  b[5] = a[7];
  b[6] = a[2];
  b[7] = a[5];
  b[8] = a[8];
}

template <typename ValueType>
void transpose4(ValueType m[16]) 
{
  swap(m[1], m[4]);
  swap(m[2], m[8]);
  swap(m[3], m[12]);
  swap(m[6], m[9]);
  swap(m[7], m[13]);
  swap(m[11], m[14]);
}

}

#endif
