#ifndef _FTK_INTERVAL_HH
#define _FTK_INTERVAL_HH

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

  static basic_interval<T> empty_interval() {
    return basic_interval<T>();
  }

  static basic_interval<T> complete_interval() {
    basic_interval<T> i;
    i.set_to_complete();
    return i;
  }

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
  bool lower_open() const {return _lower_open;}
  bool upper_open() const {return _upper_open;}

  bool complete() const {
    return lower() == -std::numeric_limits<T>::max() && 
      upper() == std::numeric_limits<T>::max();
  }

  bool empty() const {return lower() > upper();}
  bool singleton() const {return lower() == upper();}

  bool contains(T x) const {
    return x >= lower() && x <= upper();
  }

  T sample() const {
    if (singleton()) return lower();
    else return (upper() + lower()) / 2;
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
      os << i.lower() << "," << i.upper(); 
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

//////////////////////////
template <typename T>
struct disjoint_intervals {
  disjoint_intervals() {}
  disjoint_intervals(T l, T u) {
    _subintervals.insert(basic_interval<T>(l, u));
  }
  disjoint_intervals(T v) {
    _subintervals.insert(basic_interval<T>(v));
  }

  static disjoint_intervals<T> empty_interval() {
    return disjoint_intervals<T>();
  }

  static disjoint_intervals<T> complete_interval() {
    disjoint_intervals<T> I;
    I.set_to_complete();
    return I;
  }

  void set_to_empty() {
    _subintervals.clear();
  }

  void set_to_complete() {
    basic_interval<T> i;
    i.set_to_complete();
    _subintervals.clear();
    _subintervals.insert(i);
  }
 
  bool empty() const {
    return subintervals().empty();
  }

  bool singleton() const {
    return subintervals().size() == 1 && 
      subintervals().begin()->singleton();
  }

  bool contains(T x) const {
    for (const auto &I : subintervals())
      if (I.contains(x)) return true;
    return false;
  }

  bool overlaps(const basic_interval<T> &i) const {
    for (const auto &j : subintervals())
      if (j.overlaps(i)) return true;
    return false;
  }

  void join(T l, T u) {
    join(basic_interval<T>(l, u));
  }

  void join(const basic_interval<T>& i) {
    if (i.empty()) return; 

    std::list<typename std::set<basic_interval<T> >::iterator> to_erase;
    basic_interval<T> ii = i;

    for (typename std::set<basic_interval<T> >::iterator it = _subintervals.begin(); 
        it != _subintervals.end(); it ++) {
      if (it->overlaps(ii)) {
        ii.join(*it);
        to_erase.push_back(it);
      }
    }

    for (auto it : to_erase)
      _subintervals.erase(it);

    _subintervals.insert(ii);
  }
  
  void intersect(T l, T u) {
    intersect(basic_interval<T>(l, u));
  }

  void intersect(const basic_interval<T>& i) {
    disjoint_intervals<T> I;
    for (const auto &j : subintervals()) {
      auto k = i;
      k.intersect(j);
      I.join(k); 
    }
    *this = I;
  }

  friend std::ostream& operator<<(std::ostream& os, 
      const disjoint_intervals<T>& i) 
  {
    // fprintf(stderr, "%d, %d\n", i.lower(), i.upper());
    if (i.empty()) 
      os << "{empty}";
    else if (i.singleton())
      os << *i.subintervals().begin(); 
    else {
      std::stringstream ss;
      for (const auto &i : i.subintervals()) 
        ss << i << "U";
      os << ss.str().substr(0, ss.str().size()-1);
    }
    return os;
  }

  std::set<basic_interval<T> >& subintervals() {return _subintervals;}
  const std::set<basic_interval<T> >& subintervals() const {return _subintervals;}

private:
  std::set<basic_interval<T> > _subintervals;
};

}

#endif
