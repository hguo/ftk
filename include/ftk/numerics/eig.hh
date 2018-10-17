#ifndef _FTK_EIG_H
#define _FTK_EIG_H

#include <ftk/numerics/trace.hh>
#include <ftk/numerics/mulmat.hh>
#include <ftk/numerics/det.hh>
#include <ftk/numerics/cubic_solve.hh>
#include <iostream>

namespace ftk {

template <typename ValueType>
inline void eig3(ValueType m[9], ValueType eig_[3]) 
{
  using std::complex;

  complex<ValueType> eig[3];

  const ValueType b = -trace3(m),
                  c = m[4]*m[8] + m[0]*m[8] + m[0]*m[4] - m[1]*m[3] - m[5]*m[7] - m[2]*m[6],
                  d = -det3(m);

  cubic_solve(b, c, d, eig);
  for (int i=0; i<3; i++) {
    fprintf(stderr, "lambda%d=%f + i*%f\n", i, eig[i].real(), eig[i].imag());
  }
}

}

#endif
