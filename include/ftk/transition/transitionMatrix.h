#ifndef _FTK_TRANSITION_MATRIX_H
#define _FTK_TRANSITION_MATRIX_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include "ftk/transition/interval.h"
#include "ftk/transition/events.h"

class ftkTransitionMatrix {
  // friend class diy::Serialization<ftkTransitionMatrix>;
public:
  ftkTransitionMatrix() : _n0(INT_MAX), _n1(INT_MAX) {}
  // ftkTransitionMatrix(int n0, int n1) ;
  ftkTransitionMatrix(int t0, int t1, int n0, int n1) : _interval(std::make_pair(t0, t1)), _n0(n0), _n1(n1) {_matrix.resize(n0*n1);}
  ftkTransitionMatrix(ftkInterval, int n0, int n1);
  ~ftkTransitionMatrix() {}

public: // IO
  void SetToDummy() {_n0 = _n1 = 0; _matrix.clear();}
  bool Valid() const {return _matrix.size()>0;}
  inline void Print() const;
  inline void SaveAscii(const std::string& filename) const;
  
public: // modulars
  inline void Modularize();
  int NModules() const {return _lhss.size();}
  inline void GetModule(int i, std::set<int>& lhs, std::set<int>& rhs, int &event) const;

  inline void Normalize();
 
public: // access
  int& operator()(int i, int j) {return _matrix[i*n1() + j];}
  int operator()(int i, int j) const {return _matrix[i*n1() + j];}
  int& at(int i, int j) {return _matrix[i*n1() + j];}
  int at(int i, int j) const {return _matrix[i*n1() + j];}

  int t0() const {return _interval.first;} // timestep
  int t1() const {return _interval.second;}
  int n0() const {return _n0;}
  int n1() const {return _n1;}

  ftkInterval GetInterval() const {return _interval;}
  void SetInterval(const ftkInterval &i) {_interval = i;}

  inline int colsum(int j) const;
  inline int rowsum(int i) const;

private:
  // std::string MatrixFileName(const std::string& dataname, int t0, int t1) const;

private:
  ftkInterval _interval;
  int _n0, _n1;
  std::vector<int> _matrix; // match matrix

  // modulars
  std::vector<std::set<int> > _lhss, _rhss;
  std::vector<int> _events;
};



////// 
void ftkTransitionMatrix::GetModule(int i, std::set<int> &lhs, std::set<int> &rhs, int &event) const
{
  if (i>=_lhss.size()) return;

  lhs = _lhss[i];
  rhs = _rhss[i];
  event = _events[i];
}

int ftkTransitionMatrix::colsum(int j) const
{
  int sum = 0;
  for (int i=0; i<n0(); i++)
    sum += _matrix[i*n1() + j];
  return sum;
}

int ftkTransitionMatrix::rowsum(int i) const 
{
  int sum = 0;
  for (int j=0; j<n1(); j++) 
    sum += _matrix[i*n1() + j];
  return sum;
}

void ftkTransitionMatrix::Normalize()
{
  if (NModules() == 0) Modularize();
  for (int i=0; i<NModules(); i++) {
    const std::set<int> &lhs = _lhss[i];
    const std::set<int> &rhs = _rhss[i];
    for (std::set<int>::const_iterator it0=lhs.begin(); it0!=lhs.end(); it0++) 
      for (std::set<int>::const_iterator it1=rhs.begin(); it1!=rhs.end(); it1++) {
        int l = *it0, r = *it1;
        at(l, r) = 1;
      }
  }
}

void ftkTransitionMatrix::Modularize()
{
  const int n = n0() + n1();

  _lhss.clear();
  _rhss.clear();
  _events.clear();

  std::set<int> unvisited; 
  for (int i=0; i<n; i++) 
    unvisited.insert(i);

  while (!unvisited.empty()) {
    std::set<int> lhs, rhs; 
    std::vector<int> Q;
    Q.push_back(*unvisited.begin());

    while (!Q.empty()) {
      int v = Q.back();
      Q.pop_back();
      unvisited.erase(v);
      if (v<n0()) lhs.insert(v);
      else rhs.insert(v-n0());

      if (v<n0()) { // left hand side
        for (int j=0; j<n1(); j++) 
          if (at(v, j)>0 && unvisited.find(j+n0()) != unvisited.end())
            Q.push_back(j+n0());
      } else { // right hand side
        for (int i=0; i<n0(); i++) 
          if (at(i, v-n0())>0 && unvisited.find(i) != unvisited.end())
            Q.push_back(i);
      }
    }

    int event; 
    if (lhs.size() == 1 && rhs.size() == 1) {
      event = FTK_EVENT_DUMMY;
    } else if (lhs.size() == 0 && rhs.size() == 1) {
      event = FTK_EVENT_BIRTH;
    } else if (lhs.size() == 1 && rhs.size() == 0) {
      event = FTK_EVENT_DEATH;
    } else if (lhs.size() == 1 && rhs.size() == 2) {
      event = FTK_EVENT_SPLIT;
    } else if (lhs.size() == 2 && rhs.size() == 1) { 
      event = FTK_EVENT_MERGE;
    } else if (lhs.size() == 2 && rhs.size() == 2) { 
      event = FTK_EVENT_RECOMBINATION;
    } else {
      event = FTK_EVENT_COMPOUND;
    }

    _lhss.push_back(lhs);
    _rhss.push_back(rhs);
    _events.push_back(event);
  }
}

void ftkTransitionMatrix::Print() const
{
  fprintf(stderr, "Interval={%d, %d}, n0=%d, n1=%d\n", 
      _interval.first, _interval.second, _n0, _n1);

  for (int i=0; i<_n0; i++) {
    for (int j=0; j<_n1; j++) {
      fprintf(stderr, "%d\t", at(i, j));
    }
    fprintf(stderr, "\n");
  }
}

void ftkTransitionMatrix::SaveAscii(const std::string& filename) const 
{
  FILE *fp = fopen(filename.c_str(), "w");
  
  fprintf(fp, "Interval={%d, %d}, n0=%d, n1=%d\n", 
      _interval.first, _interval.second, _n0, _n1);

  for (int i=0; i<_n0; i++) {
    for (int j=0; j<_n1; j++) {
      fprintf(fp, "%d\t", at(i, j));
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}

#endif
