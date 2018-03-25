#ifndef _FEATURE_TRANSITION_MATRIX_H
#define _FEATURE_TRANSITION_MATRIX_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include "common/Interval.h"
#include "common/FeatureEvents.h"

class FeatureTransitionMatrix {
  // friend class diy::Serialization<FeatureTransitionMatrix>;
public:
  FeatureTransitionMatrix();
  FeatureTransitionMatrix(int n0, int n1);
  FeatureTransitionMatrix(int t0, int t1, int n0, int n1);
  FeatureTransitionMatrix(Interval, int n0, int n1);
  ~FeatureTransitionMatrix();

public: // IO
  void SetToDummy() {_n0 = _n1 = 0; _match.clear();}
  bool Valid() const {return _match.size()>0;}
  void Print() const;
  void SaveAscii(const std::string& filename) const;
  
public: // modulars
  void Modularize();
  int NModules() const {return _lhss.size();}
  void GetModule(int i, std::set<int>& lhs, std::set<int>& rhs, int &event) const;

  void Normalize();
 
public: // access
  int& operator()(int, int);
  int operator()(int, int) const;
  int& at(int i, int j);
  int at(int i, int j) const;

  int t0() const {return _interval.first;} // timestep
  int t1() const {return _interval.second;}
  int n0() const {return _n0;}
  int n1() const {return _n1;}

  Interval GetInterval() const {return _interval;}
  void SetInterval(const Interval &i) {_interval = i;}

  int colsum(int j) const;
  int rowsum(int i) const;

private:
  // std::string MatrixFileName(const std::string& dataname, int t0, int t1) const;

private:
  Interval _interval;
  int _n0, _n1;
  std::vector<int> _match; // match matrix

  // modulars
  std::vector<std::set<int> > _lhss, _rhss;
  std::vector<int> _events;

public:
  // vortex properties
  std::vector<float> moving_speeds; // length=n0
};


///////////
#if 0
namespace diy {
  template <> struct Serialization<FeatureTransitionMatrix> {
    static void save(diy::BinaryBuffer& bb, const FeatureTransitionMatrix& m) {
      diy::save(bb, m._interval);
      diy::save(bb, m._n0);
      diy::save(bb, m._n1);
      diy::save(bb, m._match);
      diy::save(bb, m._lhss);
      diy::save(bb, m._rhss);
      diy::save(bb, m._events);
      diy::save(bb, m.moving_speeds);
    }

    static void load(diy::BinaryBuffer&bb, FeatureTransitionMatrix& m) {
      diy::load(bb, m._interval);
      diy::load(bb, m._n0);
      diy::load(bb, m._n1);
      diy::load(bb, m._match);
      diy::load(bb, m._lhss);
      diy::load(bb, m._rhss);
      diy::load(bb, m._events);
      diy::load(bb, m.moving_speeds);
    }
  };
}
#endif

#endif
