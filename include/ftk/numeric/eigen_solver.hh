#ifndef _FTK_EIG_H
#define _FTK_EIG_H

#include <ftk/numeric/trace.hh>
#include <ftk/numeric/det.hh>
#include <ftk/numeric/matrix_multiplication.hh>
#include <ftk/numeric/cubic_solver.hh>
#include <ftk/numeric/quartic_solver.hh>
#include <ftk/numeric/vector_normalization.hh>
#include <iostream>

namespace ftk {

template <typename T>
inline void solve_eigenvalues_real_symmetric2(T m00, T m10, T m11, T eig[2])
{
  const T b = -(m00 + m11), c = m00*m11 - m10*m10;
  const T sqrt_delta = sqrt(b*b - 4*c);
  
  eig[0] = 0.5*(-b+sqrt_delta);
  eig[1] = 0.5*(-b-sqrt_delta);

  if (abs(eig[0]) < abs(eig[1]))
    std::swap(eig[0], eig[1]);
}

template <typename T>
inline void solve_eigenvalues_real_symmetric2(T m[2][2], T eig[2])
{
  return solve_eigenvalues_real_symmetric2(m[0][0], m[1][0], m[1][1], eig);
}


// compute eigenvector for a given matrix and one of its eigenvalue
template <typename T>
inline void eigvec3(const T m[9], std::complex<T> lambda, std::complex<T> v[3])
{
  std::complex<T> D = (m[0] - lambda) * (m[4] - lambda) - m[1] * m[3], 
    Dx = -m[2] * (m[4] - lambda) + m[1] * m[5],
    Dy = -m[5] * (m[0] - lambda) + m[2] * m[3];

  v[0] = Dx / D;
  v[1] = Dy / D;
  v[2] = std::complex<T>(T(1), T(0));

  vector_normalization2_3(v);
}

template <typename T>
inline void eigvec3(const T m[3][3], std::complex<T> lambda, std::complex<T> v[3])
{
  std::complex<T> D = (m[0][0] - lambda) * (m[2][2] - lambda) - m[0][1] * m[1][0], 
    Dx = -m[0][2] * (m[1][1] - lambda) + m[0][1] * m[1][2],
    Dy = -m[1][2] * (m[0][0] - lambda) + m[0][2] * m[1][0];

  v[0] = Dx / D;
  v[1] = Dy / D;
  v[2] = std::complex<T>(T(1), T(0));

  vector_normalization2_3(v);
}

template <typename T>
inline void eig2(const T m[2][2], std::complex<T> eig[2], std::complex<T> eigvec[2])
{
  const T b = -trace2(m), c = det2(m);
  qudratic_solver(b, c, eig);

  for (int i=0; i<2; i++) 
    eigvec2(m, eig[i], eigvec[i]);
}

template <typename T>
inline void eig3(const T m[9], std::complex<T> eig[3], std::complex<T> eigvec[3][3]) 
{
  // using std::complex;
  const T b = -trace3(m),
                  c = m[4]*m[8] + m[0]*m[8] + m[0]*m[4] - m[1]*m[3] - m[5]*m[7] - m[2]*m[6],
                  d = -det3(m);

  cubic_solve(b, c, d, eig);

  for (int i=0; i<3; i++) 
    eigvec3(m, eig[i], eigvec[i]);
}

template <typename T>
inline void eig3(
    const T m[3][3], 
    std::complex<T> eig[3], 
    std::complex<T> eigvec[3][3]) // untested
{
  // using std::complex;
  const T b = -trace3(m),
                  c = m[1][1]*m[2][2] + m[0][0]*m[2][2] + m[0][0]*m[1][1] - m[0][1]*m[1][0] - m[1][2]*m[2][1] - m[0][2]*m[2][0],
                  d = -det3(m);

  cubic_solve(b, c, d, eig);

  for (int i=0; i<3; i++) 
    eigvec3(m, eig[i], eigvec[i]);
}

// compute eigenvector for a given matrix and one of its eigenvalue
template <typename T>
inline void eigvec4(const T m[16], std::complex<T> lambda, std::complex<T> v[4])
{
  std::complex<T> D = det3<std::complex<T> >(m[0]-lambda, m[1], m[2], m[4], m[5]-lambda, m[6], m[8], m[9], m[10]-lambda), 
    Dx = det3<std::complex<T> >(-m[3], m[1], m[2], -m[7], m[5]-lambda, m[6], -m[11], m[9], m[10]-lambda), 
    Dy = det3<std::complex<T> >(m[0]-lambda, -m[3], m[2], m[4], -m[7], m[6], m[8], -m[11], m[10]-lambda), 
    Dz = det3<std::complex<T> >(m[0]-lambda, m[1], -m[3], m[4], m[5]-lambda, -m[7], m[8], m[9], -m[11]);

  v[0] = Dx / D;
  v[1] = Dy / D;
  v[2] = Dz / D;
  v[3] = std::complex<T>(T(1), T(0));

  vector_normalization2_4(v);
}

template <typename T>
inline void eig4(const T m[16], std::complex<T> eig[4], std::complex<T> eigvec[4][4]) 
{
  const T a0 = det4(m), 
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

  for (int i=0; i<4; i++)
    eigvec4(m, eig[i], eigvec[i]);
}

}

#endif
