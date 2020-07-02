#ifndef _FTK_ADJUGATE_HH
#define _FTK_ADJUGATE_HH

#include <ftk/ftk_config.hh>

namespace ftk {

template <typename T>
__device__ __host__
inline void adjugate2(const T M[2][2], T A[2][2])
{
  A[0][0] =  M[1][1];
  A[0][1] = -M[0][1];
  A[1][0] = -M[1][0];
  A[1][1] =  M[0][0];
}

template <typename T>
__device__ __host__
inline void adjugate3(const T M[3][3], T A[3][3])
{
  A[0][0] =  det2(M[1][1], M[1][2], M[2][1], M[2][2]);
  A[1][0] = -det2(M[1][0], M[1][2], M[2][0], M[2][2]);
  A[2][0] =  det2(M[1][0], M[1][1], M[2][0], M[2][1]);
  A[0][1] = -det2(M[0][1], M[0][2], M[2][1], M[2][2]);
  A[1][1] =  det2(M[0][0], M[0][2], M[2][0], M[2][2]);
  A[2][1] = -det2(M[0][0], M[0][1], M[2][0], M[2][1]);
  A[0][2] =  det2(M[0][1], M[0][2], M[1][1], M[1][2]);
  A[1][2] = -det2(M[0][0], M[0][2], M[1][0], M[1][2]);
  A[2][2] =  det2(M[0][0], M[0][1], M[1][0], M[1][1]);
}

}

#endif
