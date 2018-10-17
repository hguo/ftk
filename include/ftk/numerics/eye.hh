#ifndef _FTK_EYE_H
#define _FTK_EYE_H

namespace ftk {

template <typename ValueType>
void eye3(ValueType m[9]) 
{
  m[0] = 1; m[1] = 0; m[2] = 0;
  m[3] = 0; m[4] = 1; m[5] = 0;
  m[6] = 0; m[7] = 0; m[8] = 1;
}

}

#endif
