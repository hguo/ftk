#ifndef _FTK_ISNAN_H
#define _FTK_ISNAN_H

namespace ftk {

template <typename ValueType>
inline bool isnan_mat3x3(const ValueType m[])
{
  for (int i=0; i<9; i++) 
    if (isnan(m[i])) return true;
  return false;
}

template <typename ValueType>
inline bool isnan_mat3x3(const ValueType m[3][3])
{
  for (int i=0; i<3; i++) 
    for (int j=0; j<3; j++) 
      if (isnan(m[i][j])) return true;
  return false;
}

}

#endif
