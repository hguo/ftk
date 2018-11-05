#ifndef _FTK_QUARTIC_SOLVE_H
#define _FTK_QUARTIC_SOLVE_H

#include <math.h>
#include <complex>

namespace ftk {

template <typename ValueType>
inline void quartic_solve(ValueType b, ValueType c, ValueType d, ValueType e, std::complex<ValueType> x[4])
{
#if 0
  const ValueType delta = 
      256*e*e*e - 192*b*d*e*e - 128*c*c*e*e + 144*c*d*d*e - 27*d*d*d*d
    + 144*b*b*c*e*e - 6*b*b*d*d*e - 80*b*c*c*d*e + 18*b*c*d*d*d + 16*c*c*c*c*e
    - 4*c*c*c*d*d - 27*b*b*b*b*e*e + 18*b*b*b*c*d*e - 4*b*b
#endif
}

}

#endif
