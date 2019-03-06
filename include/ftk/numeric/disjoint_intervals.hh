#ifndef _FTK_DISJOINT_INTERVALS_HH
#define _FTK_DISJOINT_INTERVALS_HH

#include <random>
#include <algorithm>
#include <ftk/numeric/basic_interval.hh>

namespace ftk {

template <typename T>
struct disjoint_intervals {
  disjoint_intervals() {}
  disjoint_intervals(T l, T u) {
    _subintervals.insert(basic_interval<T>(l, u));
  }
  disjoint_intervals(T v) {
    _subintervals.insert(basic_interval<T>(v));
  }
  // disjoint_intervals(const disjoint_intervals<T>&) {};
  // disjoint_intervals(disjoint_intervals<T>&&) = default;

  // conversion from quantized intervals
  disjoint_intervals(const disjoint_intervals<long long>& I, 
      const long long factor) //  = 1000000000L) 
  {
    const T epsilon = T(1) / T(factor);
    for (auto i : I.subintervals()) { // TODO: consider open intervals
      const auto lb = i.lower_bounded() ? (i.lower() * epsilon) : (basic_interval<T>::lower_inf()), 
                 ub = i.upper_bounded() ? (i.upper() * epsilon) : (basic_interval<T>::upper_inf());
      basic_interval<T> ii(lb, ub);
      if (i.lower_open()) ii.set_lower_open();
      if (i.upper_open()) ii.set_upper_open();
      join(ii); // basic_interval<T>(lb, ub));
    }
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

  bool complete() const {
    return subintervals().size() == 1 
      && subintervals().begin()->complete();
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

  void intersect(const disjoint_intervals<T>& J) {
    disjoint_intervals<T> I;

    for (const auto &i : subintervals()) 
      for (const auto &j : J.subintervals()) {
        auto ii = i; 
        ii.intersect(j);
        // std::cerr << "ij: " << ii << std::endl;
        I.join(ii); // intersect(i, j));
      }
    *this = I;
  }

  // disjoint_intervals<T> complement() const { // TODO
  // }

  T sample() const {
    if (empty()) return std::nan("0");

    std::random_device rd;
    std::mt19937 g(rd());

    std::vector<basic_interval<T>> my_subintervals;
    for (const auto &ii : _subintervals)
      my_subintervals.push_back(ii);
    std::shuffle(my_subintervals.begin(), my_subintervals.end(), g);
    return my_subintervals[0].sample();
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
      for (const auto &ii : i.subintervals()) 
        ss << ii << "U";
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
