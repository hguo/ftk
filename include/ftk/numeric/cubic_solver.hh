#ifndef _FTK_CUBIC_SOLVER_HH
#define _FTK_CUBIC_SOLVER_HH

#include <ftk/numeric/quadratic_solver.hh>
#include <cmath>
#include <complex>

namespace ftk {

// This function returns the discriminant.  
// If disc > 0, we have one real and two complex conjugate roots.
// If disc <= 0, we have three roots; the roots may be multiple.
template <typename T>
T solve_cubic(T b, T c, T d, std::complex<T> x[3]) 
{
  T disc, q, r, dum1, s, t, term1, r13;

  q = (3.0*c - (b*b))/9.0; 
  r = (-(27.0*d) + b*(9.0*c - 2.0*(b*b)))/54.0;
  disc = q*q*q + r*r;

  x[0].imag(0); // the first root is always real.
  term1 = (b/3.0);

  if(disc > 0) { // one root real, two are complex
    s = r + sqrt(disc);
    s = ((s < 0) ? -pow(-s, (1.0/3.0)) : pow(s, (1.0/3.0)));
    t = r - sqrt(disc);
    t = ((t < 0) ? -pow(-t, (1.0/3.0)) : pow(t, (1.0/3.0)));
    x[0].real(-term1 + s + t);
    term1 += (s + t)/2.0;
    x[1].real(-term1);
    x[2].real(-term1);
    term1 = sqrt(3.0)*(-t + s)/2;
    x[1].imag(term1);
    x[2].imag(-term1);
    // return;
  } else { // the remaining options are all real
    x[1].imag(0);
    x[2].imag(0);
    if (disc == 0) { // all roots real, at least two are equal.
      r13 = ((r < 0) ? -pow(-r,(1.0/3.0)) : pow(r,(1.0/3.0)));
      x[0].real(-term1 + 2.0*r13);
      x[1].real(-(r13 + term1));
      x[2].real(-(r13 + term1));
      // return;
    } else { // all roots are real and unequal
      q = -q;
      dum1 = q*q*q;
      dum1 = acos(r/sqrt(dum1));
      r13 = 2.0*sqrt(q);
      x[0].real(-term1 + r13*cos(dum1/3.0));
      x[1].real(-term1 + r13*cos((dum1 + 2.0*M_PI)/3.0));
      x[2].real(-term1 + r13*cos((dum1 + 4.0*M_PI)/3.0));
      // return;
    }
  }
  
  return disc; 
}

template <typename T>
int solve_cubic_real(T b, T c, T d, T x[3], const T epsilon = std::numeric_limits<T>::epsilon())
{
  T disc, q, r, dum1, s, t, term1, r13;

  q = (3.0*c - (b*b))/9.0; 
  r = (-(27.0*d) + b*(9.0*c - 2.0*(b*b)))/54.0;
  disc = q*q*q + r*r;

  term1 = b / 3.0;

  if(disc > 0) { // one root real, two are complex
    s = r + sqrt(disc);
    s = ((s < 0) ? -pow(-s, (1.0/3.0)) : pow(s, (1.0/3.0)));
    t = r - sqrt(disc);
    t = ((t < 0) ? -pow(-t, (1.0/3.0)) : pow(t, (1.0/3.0)));

    x[0] = -term1 + s + t;
    return 1;
  } else { // the remaining options are all real
    if (disc == 0) { // all roots real, at least two are equal.
      r13 = ((r < 0) ? -pow(-r,(1.0/3.0)) : pow(r,(1.0/3.0)));
      x[0] = -term1 + 2.0 * r13;
      x[1] = -(r13 + term1);

      if (x[0] == x[1]) return 1;
      else return 2;
    } else { // all roots are real and unequal
      q = -q;
      dum1 = q*q*q;
      dum1 = acos(r/(sqrt(dum1)));
      r13 = 2.0*sqrt(q);
      x[0] = -term1 + r13*cos(dum1/3.0);
      x[1] = -term1 + r13*cos((dum1 + 2.0*M_PI)/3.0);
      x[2] = -term1 + r13*cos((dum1 + 4.0*M_PI)/3.0);
      return 3;
    }
  }
}


template <typename T>
void solve_cubic(const T coef[4], std::complex<T> x[3])
{
  solve_cubic(
      coef[2] / coef[3], 
      coef[1] / coef[3], 
      coef[0] / coef[3], 
      x);
}

template <typename T>
int solve_cubic_real(const T coef[4], T x[3], const T epsilon = std::numeric_limits<T>::epsilon())
{
  if (std::abs(coef[3]) <= epsilon || std::isnan(coef[3]) || std::isinf(coef[3])) {
    // fprintf(stderr, "downgrading, eps=%.20f\n", epsilon);
    return solve_quadratic_real(coef, x, epsilon);
  } else if (std::abs(coef[0]) <= epsilon || std::isnan(coef[0]) || std::isinf(coef[0])) {
    auto n = solve_quadratic_real(coef+1, x, epsilon);
    x[n++] = T(0);
    return n;
  } else {
    return solve_cubic_real(
      coef[2] / coef[3], 
      coef[1] / coef[3], 
      coef[0] / coef[3], 
      x);
  }
}

////////////////

template <typename T>
std::map<T, int> solve_cubic_real_multiplicity(T b, T c, T d)
{
  std::map<T, int> roots;
  T disc, q, r, dum1, s, t, term1, r13;

  q = (T(3)*c - (b*b))/T(9); 
  r = (-(T(27)*d) + b*(T(9)*c - T(2)*(b*b)))/T(54);
  disc = q*q*q + r*r;

  term1 = b / T(3);

  if(disc > 0) { // one root real, two are complex
    s = r + sqrt(disc);
    s = ((s < 0) ? -pow(-s, (T(1)/T(3))) : pow(s, (T(1)/T(3))));
    t = r - sqrt(disc);
    t = ((t < 0) ? -pow(-t, (T(1)/T(3))) : pow(t, (T(1)/T(3))));

    roots[-term1 + s + t] ++;
  } else { // the remaining options are all real
    if (disc == 0) { // all roots real, at least two are equal.
      r13 = ((r < 0) ? -pow(-r, (T(1)/T(3))) : pow(r, T(1)/T(3)));
      roots[-term1 + T(2) * r13] ++;
      roots[-(r13 + term1)] += 2; // multiple roots
    } else { // all roots are real and unequal
      q = -q;
      dum1 = q * q * q;
      dum1 = acos(r/(sqrt(dum1)));
      r13 = T(2) * sqrt(q);
      roots[-term1 + r13*cos(dum1/T(3))] ++;
      roots[-term1 + r13*cos((dum1 + T(2)*M_PI)/T(3))] ++;
      roots[-term1 + r13*cos((dum1 + T(4)*M_PI)/T(3))] ++;
    }
  }
  return roots;
}

template <typename T>
std::map<T, int> solve_cubic_real_multiplicity(const T P[4], const T epsilon = std::numeric_limits<T>::epsilon()) // 1e-9)
{
  if (std::abs(P[3]) <= epsilon || std::isnan(P[3]) || std::isinf(P[3])) {
    return solve_quadratic_real_multiplicity(P, epsilon);
  } else if (std::abs(P[0]) <= epsilon || std::isnan(P[0]) || std::isinf(P[0])) {
    auto roots = solve_quadratic_real_multiplicity(P+1, epsilon);
    roots[T(0)] ++;
    return roots;
  } else {
    return solve_cubic_real_multiplicity(
      P[2] / P[3], 
      P[1] / P[3], 
      P[0] / P[3]);
  }
}

}

#endif
