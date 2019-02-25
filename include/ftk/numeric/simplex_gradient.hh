#ifndef _FTK_SIMPLEX_GRADIENT_HH
#define _FTK_SIMPLEX_GRADIENT_HH

#include <ftk/numeric/matrix_inverse.hh>
#include <ftk/numeric/print.hh>

namespace ftk {

template <typename T>
void simplex2_2d_gradient1(const T X[3][2], const T f[3], T gradf[2])
{
  const T df_dmu0 = f[0] - f[2], 
          df_dmu1 = f[1] - f[2];

  const T M[2][2] = {
    {X[0][0] - X[2][0], X[1][0] - X[2][0]}, 
    {X[0][1] - X[2][1], X[1][1] - X[2][1]}
  };

  T invM[2][2];
  ftk::matrix_inverse2(M, invM);

  const T dmu0_dx0 = invM[0][0], 
          dmu1_dx1 = invM[1][0];

  gradf[0] = df_dmu0 * dmu0_dx0;
  gradf[1] = df_dmu1 * dmu1_dx1;
}

template <typename T>
void simplex2_2d_gradient2(const T X[3][2], const T f[3][2], T gradf[2][2])
{
  const T df_dmu0[2] = {f[0][0] - f[2][0], f[0][1] - f[2][1]},
          df_dmu1[2] = {f[1][0] - f[2][0], f[1][1] - f[2][1]};

  const T M[2][2] = {
    {X[0][0] - X[2][0], X[1][0] - X[2][0]}, 
    {X[0][1] - X[2][1], X[1][1] - X[2][1]}
  };

  T invM[2][2];
  const T det = ftk::matrix_inverse2(M, invM);

  const T dmu0_dx0 = invM[0][0], 
          dmu1_dx1 = invM[1][0];

  gradf[0][0] = df_dmu0[0] * dmu0_dx0;
  gradf[0][1] = df_dmu0[1] * dmu0_dx0;
  gradf[1][0] = df_dmu1[0] * dmu1_dx1;
  gradf[1][1] = df_dmu1[1] * dmu1_dx1;
}

}

#endif
