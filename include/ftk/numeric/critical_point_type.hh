#ifndef _FTK_CRITICAL_POINT_HH
#define _FTK_CRITICAL_POINT_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/eigen_solver2.hh>
#include <ftk/numeric/eigen_solver3.hh>

namespace ftk {

enum {
  CRITICAL_POINT_2D_UNKNOWN = 0,
  CRITICAL_POINT_2D_DEGENERATE = 1,
  CRITICAL_POINT_2D_ATTRACTING = 8, // 2,
  CRITICAL_POINT_2D_REPELLING = 2, // 4,
  CRITICAL_POINT_2D_SADDLE = 4, // 8,
  CRITICAL_POINT_2D_ATTRACTING_FOCUS = 16,
  CRITICAL_POINT_2D_REPELLING_FOCUS = 32,
  CRITICAL_POINT_2D_CENTER = 64,
  // for scalar field 
  CRITICAL_POINT_2D_MAXIMUM = 8, // 2,
  CRITICAL_POINT_2D_MINIMUM = 2, // 4
};

enum {
  CRITICAL_POINT_3D_UNKNOWN = 0,
  CRITICAL_POINT_3D_DEGENERATE = 1,
  CRITICAL_POINT_3D_ATTRACTING = 8,
  CRITICAL_POINT_3D_REPELLING = 2,
  CRITICAL_POINT_3D_SADDLE = 4,
  CRITICAL_POINT_3D_UNSTABLE_ATTRACTING = 16,
  CRITICAL_POINT_3D_UNSTABLE_REPELLING = 32,
  CRITICAL_POINT_3D_UNSTABLE_SADDLE = 64,
  // for scalar field
  CRITICAL_POINT_3D_MAXIMUM = 8,
  CRITICAL_POINT_3D_MINIMUM = 2,
};

template <typename T>
__host__ __device__
unsigned int critical_point_type_2d(T J[2][2], bool symmetric)
{
  if (symmetric) { // treat jacobian matrix as symmetric
    double eig[2];
    solve_eigenvalues_symmetric2x2(J, eig);
    
    if (eig[0] > 0 && eig[1] > 0) return CRITICAL_POINT_2D_MINIMUM;
    else if (eig[0] < 0 && eig[1] < 0) return CRITICAL_POINT_2D_MAXIMUM;
    else if (eig[0] * eig[1] < 0) return CRITICAL_POINT_2D_SADDLE;
    else return CRITICAL_POINT_2D_DEGENERATE;
  } else {
    std::complex<T> eig[2];
    T delta = ftk::solve_eigenvalues2x2(J, eig);
    
    if (delta >= 0) { // two real roots
      if (eig[0].real() * eig[1].real() < 0) 
        return CRITICAL_POINT_2D_SADDLE;
      else if (eig[0].real() > 0 && eig[1].real() > 0)
        return CRITICAL_POINT_2D_REPELLING;
      else if (eig[0].real() < 0 && eig[1].real() < 0)
        return CRITICAL_POINT_2D_ATTRACTING;
      else 
        return CRITICAL_POINT_2D_DEGENERATE;
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

template <typename T>
__host__ __device__
unsigned int critical_point_type_3d(T J[3][3], bool symmetric)
{
  if (symmetric) {
    T eig[3];
    solve_eigenvalues_symmetric3x3(J, eig);

    if (eig[0] * eig[1] * eig[2] == T(0)) return CRITICAL_POINT_3D_DEGENERATE;
    if (eig[0] < 0 && eig[1] < 0 && eig[2] < 0) return CRITICAL_POINT_3D_MAXIMUM;
    else if (eig[0] > 0 && eig[1] > 0 && eig[2] > 0) return CRITICAL_POINT_3D_MINIMUM;
    else return CRITICAL_POINT_2D_SADDLE;
  } else {
    // std::complex<T> eig[3];
    // ftk::solve_eigenvalues3x3(J, eig); 

    // TODO
    return CRITICAL_POINT_3D_UNKNOWN;
  }
}

inline std::string critical_point_type_to_string(const int n, const int i, bool scalar=true)
{
  if (n == 2) {
    if (i == CRITICAL_POINT_2D_UNKNOWN) return "unknown";
    else if (i == CRITICAL_POINT_2D_DEGENERATE) return "degenerate";
    else if (i == CRITICAL_POINT_2D_ATTRACTING) return scalar ? "max" : "attracting";
    else if (i == CRITICAL_POINT_2D_REPELLING) return scalar ? "min" : "repelling";
    else if (i == CRITICAL_POINT_2D_SADDLE) return "saddle";
    else if (i == CRITICAL_POINT_2D_ATTRACTING_FOCUS) return "attracting_focus";
    else if (i == CRITICAL_POINT_2D_REPELLING_FOCUS) return "repelling_focus";
    else if (i == CRITICAL_POINT_2D_CENTER) return "center";
    else return ""; // TODO: throw error
  } else if (n == 3) {
    if (i == CRITICAL_POINT_3D_UNKNOWN) return "unknown";
    else if (i == CRITICAL_POINT_3D_DEGENERATE) return "degenerate";
    else if (i == CRITICAL_POINT_3D_ATTRACTING) return scalar ? "max" : "attracting";
    else if (i == CRITICAL_POINT_3D_REPELLING) return scalar ? "min" : "repelling";
    else if (i == CRITICAL_POINT_3D_SADDLE) return "saddle";
    else if (i == CRITICAL_POINT_3D_UNSTABLE_ATTRACTING) return "unstable_attracting";
    else if (i == CRITICAL_POINT_3D_UNSTABLE_REPELLING) return "unstable_repelling";
    else if (i == CRITICAL_POINT_3D_UNSTABLE_SADDLE) return "unstable_saddle";
    else return ""; 
  } else return "";
}

}

#endif
