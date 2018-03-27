#ifndef _FTK_TRANSITION_MATRIX_H
#define _FTK_TRANSITION_MATRIX_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include "ftk/ftkInterval.h"
#include "ftk/ftkEvents.h"

class ftkTransitionMatrix {
  // friend class diy::Serialization<ftkTransitionMatrix>;
public:
  ftkTransitionMatrix();
  ftkTransitionMatrix(int n0, int n1);
  ftkTransitionMatrix(int t0, int t1, int n0, int n1);
  ftkTransitionMatrix(ftkInterval, int n0, int n1);
  ~ftkTransitionMatrix();

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

  ftkInterval GetInterval() const {return _interval;}
  void SetInterval(const ftkInterval &i) {_interval = i;}

  int colsum(int j) const;
  int rowsum(int i) const;

private:
  // std::string MatrixFileName(const std::string& dataname, int t0, int t1) const;

private:
  ftkInterval _interval;
  int _n0, _n1;
  std::vector<int> _match; // match matrix

  // modulars
  std::vector<std::set<int> > _lhss, _rhss;
  std::vector<int> _events;

public:
  // vortex properties
  // std::vector<float> moving_speeds; // length=n0
};

#endif
