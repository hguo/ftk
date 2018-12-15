#include <ftk/numeric/polynomial.hh>
#include <ftk/numeric/print.hh>
#include <iostream>

int main(int argc, char **argv)
{
  float p0[3] = {1, 1, 1}, 
        p1[2] = {1, 1};
  float p[4];

  ftk::print_polynomial("p0", p0, 2);
  ftk::print_polynomial("p1", p1, 1);

  // ftk::polynomial_multiplication(p0, 2, p1, 1, p);
  ftk::polynomial_add_in_place(p0, 2, p1, 1);
  ftk::polynomial_scalar_multiplication(p0, 2.f, 2.f);
  ftk::print_polynomial("p0", p0, 2);

  fprintf(stderr, "eval_p0(2)=%f\n", ftk::polynomial_evaluate(p0, 2, 2.f));

  ftk::polynomial_derivative(p0, 2, p);
  ftk::print_polynomial("p", p, 1);

  fprintf(stderr, "eval_p(2)=%f\n", ftk::polynomial_derivative_evaluate(p0, 2, 2.f));
  fprintf(stderr, "eval_p(2)=%f\n", ftk::polynomial_evaluate(p, 1, 2.f));

  return 0;
}
