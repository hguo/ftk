#ifndef _FTK_PARALLEL_VECTOR_H
#define _FTK_PARALLEL_VECTOR_H

#include <ftk/numerics/invmat.hh>
#include <ftk/numerics/eig.hh>
#include <ftk/numerics/transpose.hh>
#include <ftk/numerics/isnan.hh>

namespace ftk {

template <typename ValueType>
inline int solve_parallel_vector_barycentric(const ValueType V[9], const ValueType W[9], ValueType lambda[9])
{
  ValueType invV[9], invW[9];
  const auto detV = invmat3(V, invV),
             detW = invmat3(W, invW);

  if (isnan_mat3x3(invW) || isnan_mat3x3(invV)) return false;
  // if (detW < 1e-4) return false;

  ValueType m[9];
  mulmat3(invW, V, m);
  // if (detW > detV) mulmat3(invW, V, m); 
  // else mulmat3(invV, W, m); // very wrong...

  std::complex<ValueType> eig[3], eigvec[3][3];
  eig3(m, eig, eigvec);

  // fprintf(stderr, "eigs_impl0: (%f, %f), (%f, %f), (%f, %f)\n", 
  //     eig[0].real(), eig[0].imag(), eig[1].real(), eig[1].imag(), eig[2].real(), eig[2].imag());

  int n_solutions = 0;
  for (int i=0; i<3; i++) {
    if (isnan(eig[i].real()) || eig[i].imag() != 0) continue; // non-real eigenvalue
    ValueType l[3] = {eigvec[i][0].real(), eigvec[i][1].real(), eigvec[i][2].real()};
    if (l[0] < 0 || l[1] < 0) continue;
    const auto sum = l[0] + l[1] + l[2];
    lambda[n_solutions*3+0] = l[0] / sum;
    lambda[n_solutions*3+1] = l[1] / sum;
    lambda[n_solutions*3+2] = l[2] / sum;
    if (isnan(lambda[0]) || isnan(lambda[1]) || isnan(lambda[2])) continue;
    n_solutions ++;
    // return true;
  }

  if (n_solutions)
    fprintf(stderr, "n_solutions=%d\n", n_solutions);
  return n_solutions;
  // return false;
}

template <typename ValueType>
inline int solve_parallel_vector_barycentric(const ValueType V[3][3], const ValueType W[3][3], ValueType lambda[3][3])
{
  ValueType VV[] = {V[0][0], V[0][1], V[0][2], V[1][0], V[1][1], V[1][2], V[2][0], V[2][1], V[2][2]};
  ValueType WW[] = {W[0][0], W[0][1], W[0][2], W[1][0], W[1][1], W[1][2], W[2][0], W[2][1], W[2][2]};
  transpose3(VV);
  transpose3(WW);
  
  ValueType lambda1[9];
  const auto n = solve_parallel_vector_barycentric(VV, WW, lambda1);
  for (int i = 0; i < n; i ++)
    for (int j = 0; j < 3; j ++)
      lambda[i][j] = lambda1[i*3+j];

  return n;
}

template <typename ValueType>
inline bool solve_parallel_vector_bilinear(const ValueType V[12], const ValueType W[12], ValueType lambda[2])
{
  // TODO
}

}
#endif
