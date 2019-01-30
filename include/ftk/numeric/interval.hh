#ifndef _FTK_CLOSED_INTERVAL_HH
#define _FTK_CLOSED_INTERVAL_HH

#include <limits>
#include <algorithm>
#include <iostream>

namespace ftk {

template <typename T>
struct basic_interval {
  basic_interval() {set_to_empty();}
  basic_interval(T l, T u) : _lower(l), _upper(u) {}
  basic_interval(T v) : _lower(v), _upper(v) {}

  void set_to_empty() {
    _lower = std::numeric_limits<T>::max();
    _upper = -std::numeric_limits<T>::max();
  }

  void set_to_complete() {
    _lower = -std::numeric_limits<T>::max();
    _upper = std::numeric_limits<T>::max();
  }

  void set_to_singleton(T x) {
    _lower = x;
    _upper = x;
  }

  T lower() const {return _lower;}
  T upper() const {return _upper;}

  bool complete() const {
    return lower() == -std::numeric_limits<T>::max() && 
      upper() == std::numeric_limits<T>::max();
  }

  bool empty() const {return lower() > upper();}
  bool singleton() const {return lower() == upper();}

  bool contains(T x) const {
    return x >= lower() && x <= upper();
  }

  bool overlaps(const basic_interval<T> &i) const {
    if (lower() > i.upper() || upper() < i.lower()) return false;
    else return true;
  }

  bool join(const basic_interval<T>& i) {
    if (overlaps(i)) {
      _lower = std::min(lower(), i.lower());
      _upper = std::max(upper(), i.upper());
      return true;
    }
    else return false; // not possible because there is no overlaps
  }

  void intersect(const basic_interval<T> &i) {
    if (overlaps(i)) {
      _lower = std::max(lower(), i.lower());
      _upper = std::min(upper(), i.upper());
    } else {
      set_to_empty();
    }
  }

  friend basic_interval<T> join(const basic_interval<T> &i, 
      const basic_interval<T> &j)
  {
    basic_interval<T> k;
    if (i.overlaps(j)) {
      k = i;
      k.join(j);
    }
    return k;
  }

  friend basic_interval<T> intersect(const basic_interval<T> &i, 
      const basic_interval<T> &j)
  {
    basic_interval<T> k = i;
    k.intersect(j);
    return k;
  }

  friend std::ostream& operator<<(std::ostream& os, 
      const basic_interval<T>& i) 
  {
    // fprintf(stderr, "%d, %d\n", i.lower(), i.upper());
    if (i.empty()) 
      os << "{empty}";
    else if (i.singleton())
      os << "{" << i.lower() << "}";
    else 
      os << "[" << i.lower() << "," << i.upper() << "]";
    return os;
  }

private:
  T _lower, _upper;
};

}

#endif
