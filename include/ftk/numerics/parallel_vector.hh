#ifndef _FTK_PARALLEL_VECTOR_H
#define _FTK_PARALLEL_VECTOR_H

#include <ftk/numerics/invmat.hh>
#include <ftk/numerics/eig.hh>

namespace ftk {

template <typename ValueType>
inline bool parallel_vector(const ValueType V[9], const ValueType W[9], ValueType lambda[3])
{
  ValueType invV[9], invW[9];
  const auto detV = invmat3(V, invV),
             detW = invmat3(W, invW);

  if (detW < 1e-4) return false;

  ValueType m[9];
  mulmat3(invW, V, m);
  // if (detW > detV) mulmat3(invW, V, m);
  // else mulmat3(invV, W, m);

  std::complex<ValueType> eig[3], eigvec[3][3];
  eig3(m, eig, eigvec);

  for (int i=0; i<3; i++) {
    if (eig[i].imag() != 0) continue; // non-real eigenvalue
    ValueType l[3] = {eigvec[i][0].real(), eigvec[i][1].real(), eigvec[i][2].real()};
    if (l[0] < 0 || l[1] < 0) continue;
    const auto sum = l[0] + l[1] + l[2];
    lambda[0] = l[0] / sum;
    lambda[1] = l[1] / sum;
    lambda[2] = l[2] / sum;
    return true;
  }

  return false;
}

}

#endif
