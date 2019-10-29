#ifndef _FTK_CRITICAL_POINT_HH
#define _FTK_CRITICAL_POINT_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/eigen_solver2.hh>

namespace ftk {

enum {
  CRITICAL_POINT_2D_UNKNOWN = 0,
  CRITICAL_POINT_2D_DEGENERACY = 0x1,
  CRITICAL_POINT_2D_ATTRACTING = 0x10,
  CRITICAL_POINT_2D_REPELLING = 0x100,
  CRITICAL_POINT_2D_SADDLE = 0x1000,
  CRITICAL_POINT_2D_ATTRACTING_FOCUS = 0x10000,
  CRITICAL_POINT_2D_REPELLING_FOCUS = 0x100000,
  CRITICAL_POINT_2D_CENTER = 0x1000000,
  // for scalar field 
  CRITICAL_POINT_2D_MINIMUM = 0x10,
  CRITICAL_POINT_2D_MAXIMUM = 0x100,
};

template <typename T>
unsigned int critical_point_type_2d(T J[2][2], bool symmetric)
{
  if (symmetric) { // treat jacobian matrix as symmetric
    double eig[2];
    solve_eigenvalues_symmetric2x2(J, eig);
    
    if (eig[0] > 0 && eig[1] > 0) return CRITICAL_POINT_2D_MAXIMUM;
    else if (eig[0] < 0 && eig[1] < 0) return CRITICAL_POINT_2D_MINIMUM;
    else if (eig[0] * eig[1] < 0) return CRITICAL_POINT_2D_SADDLE;
    else return CRITICAL_POINT_2D_DEGENERACY;
  } else {
    std::complex<double> eig[2];
    double delta = ftk::solve_eigenvalues2x2(J, eig);
    
    if (delta >= 0) { // two real roots
      if (eig[0].real() * eig[1].real() < 0) 
        return CRITICAL_POINT_2D_SADDLE;
      else if (eig[0].real() > 0 && eig[1].real() > 0)
        return CRITICAL_POINT_2D_REPELLING;
      else if (eig[0].real() < 0 && eig[1].real() < 0)
        return CRITICAL_POINT_2D_ATTRACTING;
      else 
        return CRITICAL_POINT_2D_DEGENERACY;
    } else { // two conjugate roots
      if (eig[0].real() < 0) 
        return CRITICAL_POINT_2D_ATTRACTING_FOCUS;
      else if (eig[0].real() > 0) 
        return CRITICAL_POINT_2D_REPELLING_FOCUS;
      else 
        return CRITICAL_POINT_2D_CENTER;
    }
  }
}

}

#endif
