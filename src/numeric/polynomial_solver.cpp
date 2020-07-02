#include <iostream>
#include <cassert>
#include "ftk/numeric/polynomial_solver.hh"

#if FTK_HAVE_MPSOLVE
#include <mps/mps.h>
#include <gmp.h>
#endif

namespace ftk {

template <>
bool solve_polynomials(const double * x, int n, double * root_real, double * root_im)
{
#if FTK_HAVE_MPSOLVE
	mps_context * s = mps_context_new();
	mps_context_select_algorithm(s, MPS_ALGORITHM_SECULAR_GA);
	mps_monomial_poly * poly = mps_monomial_poly_new(s, n);
	for(int i=0; i<=n; i++){
		mps_monomial_poly_set_coefficient_d(s, poly, i, x[i], 0);
	}
	mps_context_set_input_poly(s, MPS_POLYNOMIAL(poly));
	mps_mpsolve (s);
	cplx_t *results = cplx_valloc(n);
	mps_context_get_roots_d(s, &results, NULL);
	for (int i=0; i<n; i++){
		root_real[i] = cplx_Re(results[i]);
		root_im[i] = cplx_Im(results[i]);
	}
	cplx_vfree(results);
	if(poly) mps_monomial_poly_free(s, MPS_POLYNOMIAL(poly));
	mps_context_free(s);
	return true;
#else 
  fprintf(stderr, "[ftk] FATAL: no polynomial solvers are linked.\n");
  assert(false);
  return false;
#endif
}

}
