#ifndef _FTK_EIG_H
#define _FTK_EIG_H

#include <ftk/numerics/trace.hh>
#include <ftk/numerics/det.hh>
#include <ftk/numerics/mulmat.hh>
#include <ftk/numerics/cubic_solve.hh>
#include <ftk/numerics/quartic_solve.hh>
#include <iostream>

namespace ftk {

// compute eigenvector for a given matrix and one of its eigenvalue
template <typename ValueType>
inline void eigvec3(const ValueType m[9], std::complex<ValueType> lambda, std::complex<ValueType> v[3])
{
  std::complex<ValueType> D = (m[0] - lambda) * (m[4] - lambda) - m[1] * m[3], 
    Dx = -m[2] * (m[4] - lambda) + m[1] * m[5],
    Dy = -m[5] * (m[0] - lambda) + m[2] * m[3];

  v[0] = Dx / D;
  v[1] = Dy / D;
  v[2] = std::complex<ValueType>(ValueType(1), ValueType(0));
}

template <typename ValueType>
inline void eig3(const ValueType m[9], std::complex<ValueType> eig[3], std::complex<ValueType> eigvec[3][3]) 
{
  // using std::complex;
  const ValueType b = -trace3(m),
                  c = m[4]*m[8] + m[0]*m[8] + m[0]*m[4] - m[1]*m[3] - m[5]*m[7] - m[2]*m[6],
                  d = -det3(m);

  cubic_solve(b, c, d, eig);

  for (int i=0; i<3; i++) 
    eigvec3(m, eig[i], eigvec[i]);

#if 0
  fprintf(stderr, "matrix:\n[[%f, %f, %f], [%f, %f, %f], [%f, %f, %f]]\neigenvalues:\n[%f+j%f, %f+j%f, %f+j%f]\n", 
      m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], 
      eig[0].real(), eig[0].imag(), eig[1].real(), eig[1].imag(), eig[2].real(), eig[2].imag());
  fprintf(stderr, "eigs_impl0: (%f, %f), (%f, %f), (%f, %f)\n", 
      eig[0].real(), eig[0].imag(), eig[1].real(), eig[1].imag(), eig[2].real(), eig[2].imag());
  for (int k=0; k<3; k++)
    fprintf(stderr, "eigenvector%d: %f+j%f, %f+j%f, %f+j%f\n", k,
        eigvec[k][0].real(), eigvec[k][0].imag(), eigvec[k][1].real(), eigvec[k][1].imag(), eigvec[k][2].real(), eigvec[k][2].imag());
#endif
}

template <typename ValueType>
inline void eig4(const ValueType m[16], std::complex<ValueType> eig[4], std::complex<ValueType> eigvec[4][4]) 
{
  const ValueType a0 = det4(m), 
                  a1 =  m[0]  * m[11] * m[14] - m[0]  * m[10] * m[15] - m[11] * m[12] * m[2] + m[10] * m[12] * m[3]
                      + m[1]  * m[10] * m[4]  + m[1]  * m[15] * m[4]  - m[13] * m[3]  * m[4] - m[0]  * m[10] * m[5]
                      + m[11] * m[14] * m[5]  - m[0]  * m[15] * m[5]  - m[10] * m[15] * m[5] + m[12] * m[3]  * m[5]
                      - m[11] * m[13] * m[6]  - m[1]  * m[12] * m[7]  + m[0]  * m[13] * m[7] + m[10] * m[13] * m[7]
                      + m[15] * m[2]  * m[8]  - m[14] * m[3]  * m[8]  + m[2]  * m[5]  * m[8] - m[1]  * m[6]  * m[8]
                      - m[2]  * m[4]  * m[9]  + m[0]  * m[6]  * m[9]  + m[15] * m[6]  * m[9] - m[14] * m[7]  * m[9],
                  a2 =  m[0]  * m[10] - m[11] * m[14] + m[0]  * m[15] + m[10] * m[15]
                      - m[12] * m[3]  - m[1]  * m[4]  + m[0]  * m[5]  + m[10] * m[5]
                      + m[15] * m[5]  - m[13] * m[7]  - m[2]  * m[8]  - m[6]  * m[9],
                  a3 = -trace4(m);

  quartic_solve(a3, a2, a1, a0, eig); 

  // TODO
}

}

#endif
