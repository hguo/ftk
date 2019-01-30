#ifndef _FTK_BASIC_INTERVAL_HH
#define _FTK_BASIC_INTERVAL_HH

#include <limits>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <set>
#include <list>

namespace ftk {

template <typename T>
struct basic_interval {
  basic_interval() {set_to_empty();}
  basic_interval(T l, T u) : _lower(l), _upper(u) {}
  basic_interval(T v) : _lower(v), _upper(v) {}

  static T upper_inf() {
    if (std::numeric_limits<T>::has_infinity)
      return std::numeric_limits<T>::infinity();
    else 
      return std::numeric_limits<T>::max();
  }

  static T lower_inf() {
    if (std::numeric_limits<T>::has_infinity)
      return -std::numeric_limits<T>::infinity();
    else 
      return std::numeric_limits<T>::min();
  }

  void set_lower_open() {_lower_open = true;}
  void set_lower_closed() {_lower_open = false;}
  void set_upper_open() {_upper_open = true;}
  void set_upper_closed() {_upper_open = false;}
  void set_open() {
    set_lower_open();
    set_upper_open();
  }
  void set_closed() {
    set_lower_closed();
    set_upper_closed();
  }

  static basic_interval<T> empty_interval() {
    return basic_interval<T>();
  }

  static basic_interval<T> complete_interval() {
    basic_interval<T> i;
    i.set_to_complete();
    return i;
  }

  void set_to_empty() {
    _lower = lower_inf();
    _upper = upper_inf(); 
  }

  void set_to_complete() {
    _lower = lower_inf();
    _upper = upper_inf();
  }

  void set_to_singleton(T x) {
    _lower = x;
    _upper = x;
  }

  T lower() const {return _lower;}
  T upper() const {return _upper;}
  bool lower_open() const {
    if (lower_unbounded()) return true;
    else return _lower_open;
  }
  bool upper_open() const {
    if (upper_unbounded()) return true; 
    else return _upper_open;
  }
  bool lower_bounded() const {
    return lower() != lower_inf();
  }
  bool upper_bounded() const {
    return upper() != upper_inf();
  }
  bool lower_unbounded() const {
    return !lower_bounded();
  }
  bool upper_unbounded() const {
    return !upper_bounded();
  }

  bool complete() const {
    return lower() == lower_inf() && 
      upper() == upper_inf(); 
  }

  bool empty() const {return lower() > upper();}
  bool singleton() const {return lower() == upper();}

  bool contains(T x) const {
    return x >= lower() && x <= upper();
  }

  T sample() const {
    if (complete()) return T(0);
    else if (singleton()) return lower();
    else if (upper_unbounded()) return lower() + T(1000000L);
    else if (lower_unbounded()) return upper() - T(1000000L);
    else return (upper() + lower()) / T(2);
  }

  bool overlaps(const basic_interval<T> &i) const { // FIXME: open intervals
    if (lower() > i.upper() || upper() < i.lower()) return false;
    else return true;
  }

  bool join(const basic_interval<T>& i) { // FIXME: open intervals
    if (overlaps(i)) {
      _lower = std::min(lower(), i.lower());
      _upper = std::max(upper(), i.upper());
      return true;
    }
    else return false; // not possible because there is no overlaps
  }

  void intersect(const basic_interval<T> &i) { // FIXME: open intervals
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

  friend basic_interval<T> intersect(
      const basic_interval<T> &i, 
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
    else {
      if (i.lower_open()) os << "(";
      else os << "[";

      if (i.lower_unbounded()) os << "-inf";
      else os << i.lower(); 
      os << ",";

      if (i.upper_unbounded()) os << "+inf";
      else os << i.upper();

      if (i.upper_open()) os << ")";
      else os << "]";
    }
    return os;
  }

  friend bool operator<(const basic_interval<T>& i, 
      const basic_interval<T>& j)
  {
    return i.lower() < j.lower();
  }

private:
  T _lower, _upper;
  bool _lower_open=false, _upper_open=false;
};

}

#endif
