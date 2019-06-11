#ifndef _FTK_POLYNOMIAL_SOLVER_HH
#define _FTK_POLYNOMIAL_SOLVER_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/polynomial.hh>

namespace ftk {

template <typename T>
bool solve_polynomials(const T * x, int n, double * root_real, double * root_im);

}

#endif
