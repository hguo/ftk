#ifndef _FTK_FIXED_POINT_HH
#define _FTK_FIXED_POINT_HH

namespace ftk {

template <typename I=long long/*long integer type*/, int factor=65536/*16bits*/>
struct fixed_point {
  fixed_point() : num(0) {}
  fixed_point(const fixed_point& x) : num(x.num) {}
  template <typename J> fixed_point(const J& x) {
    num = static_cast<I>(factor * x);
  }

  // explicit fixed_point(int x) {num = static_cast<I>(factor * x);}

  fixed_point& operator=(const fixed_point& x) {num = x.num; return *this;}
  template <typename J>
  fixed_point& operator=(const J& x) {
    num = static_cast<I>(factor * x);
    return *this;
  }

  I integer() const {return num;}

  fixed_point& operator+=(const fixed_point& x) {
    num += x.num;
    return *this;
  }
  fixed_point& operator-=(const fixed_point& x) {
    num -= x.num;
    return *this;
  }
  fixed_point& operator*=(const fixed_point& x) {
    num = (num * x.num) / factor;
    return *this;
  }
  fixed_point& operator/=(const fixed_point& x) {
    num = num / x.num / factor;
    return *this;
  }

  fixed_point& operator++() {num += factor; return *this;}
  fixed_point& operator--() {num -= factor; return *this;}

  bool operator!() const {return !num;}

  bool operator<(const fixed_point& x) const {return num < x.num;}
  bool operator==(const fixed_point& x) const {return num == x.num;}
  bool operator>(const fixed_point& x) const {return num > x.num;}

private:
  I num;
};

template <typename I, int factor>
fixed_point<I, factor> operator+(const fixed_point<I, factor>& a, const fixed_point<I, factor>& b)
{
  fixed_point<I> t(a);
  return t += b;
}

template <typename I, int factor>
fixed_point<I, factor> operator-(const fixed_point<I, factor>& a, const fixed_point<I, factor>& b)
{
  fixed_point<I> t(a);
  return t -= b;
}

template <typename I, int factor>
fixed_point<I, factor> operator-(const fixed_point<I, factor>& a)
{
  return static_cast<fixed_point<I, factor>>(0) - a;
}

template <typename I, int factor>
fixed_point<I, factor> operator*(const fixed_point<I, factor>& a, const fixed_point<I, factor>& b)
{
  fixed_point<I> t(a);
  return t *= b;
}

template <typename I, int factor>
std::ostream& operator<<(std::ostream& os, const fixed_point<I, factor>& x)
{
  os << x.integer() << "/" << factor;
  return os;
}

} // ftk

#endif
