#ifndef _FTK_RATIONAL_NUMBER_HH
#define _FTK_RATIONAL_NUMBER_HH

#include <ftk/ftk_config.hh>
#include <ftk/numeric/gcd.hh>
#include <iostream>
#include <iomanip>

namespace ftk {

// simplified from boost, adding __device__ __host__ flags to enable running on GPU

template <typename I=long long> // I is integer type
class rational
{
public:
  __device__ __host__ rational() : num(0), den(1) {} // zero
  __device__ __host__ rational(I n) : num(n), den(1) {} // n/1
  __device__ __host__ rational(I n, I d) : num(n), den(d) {} // n/d
  template <typename J> __device__ __host__
  rational(const rational<J>& r);

  __device__ __host__ rational& operator=(I n) {num = n; den = 1;}
  __device__ __host__ rational& assign(I n, I d) {num = n; den = d;}

  __device__ __host__ I numerator() const {return num;}
  __device__ __host__ I denominator() const {return den;}

  // arithmetic operators
  __device__ __host__ rational& operator+= (const rational& r) {
    I g = gcd(den, r.den);
    den /= g;  // = b1 from the calculations above
    num = num * (r.den / g) + r.num * den;
    g = gcd(num, g);
    num /= g;
    den *= r.den/g;
    return *this;
  }
  __device__ __host__ rational& operator-= (const rational& r) {
    I g = gcd(den, r.den);
    den /= g;
    num = num * (r.den / g) - r.num * den;
    g = gcd(num, g);
    num /= g;
    den *= r.den/g;
    return *this;
  }
  __device__ __host__ rational& operator*= (const rational& r) {
    I g1 = gcd(num, r.den);
    I g2 = gcd(r.num, den);
    num = (num/g1) * (r.num/g2);
    den = (den/g2) * (r.den/g1);
    return *this;
  }
  __device__ __host__ rational& operator/= (const rational& r) {
  	// assert(r.num != T(0))
    if (num == I(0)) return *this;

    I gcd1 = gcd(num, r.num);
    I gcd2 = gcd(r.den, den);
    num = (num/gcd1) * (r.den/gcd2);
    den = (den/gcd2) * (r.num/gcd1);

    if (den < I(0)) {
        num = -num;
        den = -den;
    }
    return *this;
  }
  
  // arithmetic/integer operators
  __device__ __host__ rational& operator+= (I i) {num += i * den; return *this;}
  __device__ __host__ rational& operator-= (I i) {num -= i * den; return *this;}
  __device__ __host__ rational& operator*= (I i) {I g = gcd(i, den); num *= i / g; den /= g; return *this;}
  __device__ __host__ rational& operator/= (I i) {
  	// assert(i != I(0));
  	if (num == I(0)) return *this;
  	I g = gcd(num, i);
  	num /= g;
  	den *= i / g;
  	if (den < I(0)) {
  	  num = -num;
  	  den = -den;
  	}
  	return *this;
  }

  // inc/dec
  __device__ __host__ const rational& operator++() {num += den; return *this;}
  __device__ __host__ const rational& operator--() {num -= den; return *this;}

  // not
  __device__ __host__ bool operator!() const {return !num;}

  // booltype
  // __device__ __host__ operator bool_type() const;

  // comparison
  __device__ __host__ bool operator< (const rational& r) const {
    struct { I  n, d, q, r; }
			ts = { this->num, this->den, static_cast<I>(this->num / this->den),
			static_cast<I>(this->num % this->den) },
			rs = { r.num, r.den, static_cast<I>(r.num / r.den),
			static_cast<I>(r.num % r.den) };
    unsigned  reverse = 0u;

    while ( ts.r < I(0) )  { ts.r += ts.d; --ts.q; }
    while ( rs.r < I(0) )  { rs.r += rs.d; --rs.q; }

    for ( ;; )
    {
      if ( ts.q != rs.q ) return reverse ? ts.q > rs.q : ts.q < rs.q;
      reverse ^= 1u;

      if ( (ts.r == I(0)) || (rs.r == I(0)) ) break;

      ts.n = ts.d;         ts.d = ts.r;
      ts.q = ts.n / ts.d;  ts.r = ts.n % ts.d;
      rs.n = rs.d;         rs.d = rs.r;
      rs.q = rs.n / rs.d;  rs.r = rs.n % rs.d;
    }

    if ( ts.r == rs.r ) return false;
    else return ( ts.r != I(0) ) != static_cast<bool>( reverse );
  }
  __device__ __host__ bool operator== (const rational& r) const {
    return ((num == r.num) && (den == r.den));
  }

  // comparison w/ integers
  __device__ __host__ bool operator< (I i) const {
  	I q = num / den, r = num % den;
  	while (r < I(0)) {r += den; --q;}
  	return q < i;
  }
  __device__ __host__ bool operator> (I i) const {
  	return operator==(i) ? false : !operator<(i);
  }
  __device__ __host__ bool operator== (I i) const {
  	return ((den == I(1)) && (num == i));
  }

private:
  __device__ __host__ void normalize() {
    // assert(den != I(0));
    if (num == I(0)) {
      den = I(1);
      return;
    }

    I g = gcd(num, den);
    num /= g;
    den /= g;

    if (den < I(0)) {
      num = -num;
      den = -den;
    }
  }

private:
  I num, den; // numerator and denominator;
};

// Unary operators
template <typename I> __device__ __host__ rational<I> 
operator+ (const rational<I>& a, const rational<I> &b) {
  rational<I> t(b);
  return t += b;
}
template <typename I> __device__ __host__ rational<I> 
operator- (const rational<I>& a, const rational<I> &b) {
  rational<I> t(b);
  return t -= b;
}

// Reversed order operators for - and / between (types convertible to) I and rational
template <typename I, typename II> __device__ __host__ 
inline rational<I> operator- (II i, const rational<I>& r); 
template <typename I, typename II> __device__ __host__ 
inline rational<I> operator/ (II i, const rational<I>& r); 

// Absolute value
template <typename I> __device__ __host__
rational<I> abs (const rational<I>& r) {
  return r.numerator() >= I(0) ? r: -r;
}

// Input and output
template <typename I> std::istream& operator>> (std::istream& is, rational<I>& r);
template <typename I> std::ostream& operator<< (std::ostream& os, const rational<I>& r) {
  std::ostringstream ss;

  ss.copyfmt( os );
  ss.tie( NULL );
  ss.exceptions( std::ios::goodbit );
  ss.width( 0 );
  ss << std::noshowpos << std::noshowbase << '/' << r.denominator();

  // The numerator holds the showpos, internal, and showbase flags.
  std::string const   tail = ss.str();
  std::streamsize const  w =
      os.width() - static_cast<std::streamsize>( tail.size() );

  ss.clear();
  ss.str( "" );
  ss.flags( os.flags() );
  ss << std::setw( w < 0 || (os.flags() & std::ios::adjustfield) !=
                   std::ios::internal ? 0 : w ) << r.numerator();
  return os << ss.str() + tail;
}

// Type conversion
template <typename T, typename I> __device__ __host__
T rational_cast (const rational<I>& r){
  return static_cast<T>(r.numerator())/static_cast<T>(r.denominator());
}

}// namespace ftk

#endif
