#ifndef _FTK_FIXED_POINT_HH
#define _FTK_FIXED_POINT_HH

#include <ftk/ftk_config.hh>

namespace ftk {

template <typename I=long long/*long integer type*/, int factor=FTK_FP_PRECISION>
struct fixed_point {
  __device__ __host__ fixed_point() : num(0) {}
  __device__ __host__ fixed_point(const fixed_point& x) : num(x.num) {}
  template <typename J> __device__ __host__ fixed_point(const J& x) {
    num = static_cast<I>(factor * x);
  }

  fixed_point& __device__ __host__ operator=(const fixed_point& x) {num = x.num; return *this;}
  template <typename J>
  fixed_point& __device__ __host__ operator=(const J& x) {
    num = static_cast<I>(factor * x);
    return *this;
  }

  I __device__ __host__ integer() const {return num;}
  int __device__ __host__ to_int() const {return static_cast<int>(num / factor);}
  float __device__ __host__ to_float() const {return static_cast<float>(num) / factor;}
  double __device__ __host__ to_double() const {return static_cast<double>(num) / factor;}

  fixed_point& __device__ __host__ operator+=(const fixed_point& x) {
    num += x.num;
    return *this;
  }
  fixed_point& __device__ __host__ operator-=(const fixed_point& x) {
    num -= x.num;
    return *this;
  }
  fixed_point& __device__ __host__ operator*=(const fixed_point& x) {
    num = (num * x.num) / factor;
    return *this;
  }
  fixed_point& __device__ __host__ operator/=(const fixed_point& x) {
    num = num / x.num / factor;
    return *this;
  }

  fixed_point& __device__ __host__ operator++() {num += factor; return *this;}
  fixed_point& __device__ __host__ operator--() {num -= factor; return *this;}

  bool __device__ __host__ operator!() const {return !num;}

  bool __device__ __host__ operator<(const fixed_point& x) const {return num < x.num;}
  bool __device__ __host__ operator<=(const fixed_point& x) const {return num <= x.num;}
  bool __device__ __host__ operator==(const fixed_point& x) const {return num == x.num;}
  bool __device__ __host__ operator>(const fixed_point& x) const {return num > x.num;}
  bool __device__ __host__ operator>=(const fixed_point& x) const {return num >= x.num;}

private:
  I num;
};

template <typename I, int factor>
fixed_point<I, factor> __device__ __host__ operator+(const fixed_point<I, factor>& a, const fixed_point<I, factor>& b)
{
  fixed_point<I> t(a);
  return t += b;
}

template <typename I, int factor>
fixed_point<I, factor> __device__ __host__ operator-(const fixed_point<I, factor>& a, const fixed_point<I, factor>& b)
{
  fixed_point<I> t(a);
  return t -= b;
}

template <typename I, int factor>
fixed_point<I, factor> __device__ __host__ operator-(const fixed_point<I, factor>& a)
{
  return static_cast<fixed_point<I, factor>>(0) - a;
}

template <typename I, int factor>
fixed_point<I, factor> __device__ __host__ operator*(const fixed_point<I, factor>& a, const fixed_point<I, factor>& b)
{
  fixed_point<I> t(a);
  return t *= b;
}

template <typename I, int factor>
fixed_point<I, factor> __device__ __host__ operator/(const fixed_point<I, factor>& a, const fixed_point<I, factor>& b)
{
  fixed_point<I> t(a);
  return t /= b;
}

template <typename I, int factor>
std::ostream& __device__ __host__ operator<<(std::ostream& os, const fixed_point<I, factor>& x)
{
  os << x.integer() << "/" << factor;
  return os;
}

} // ftk

#endif
