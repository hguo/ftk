#ifndef _FTK_TRANSITION_MATRIX_H
#define _FTK_TRANSITION_MATRIX_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include "ftk/transition/interval.h"
#include "ftk/transition/event.h"

class ftkTransitionMatrix {
public:
  ftkTransitionMatrix() {};
  ftkTransitionMatrix(int t0, int t1, int n0, int n1) : _interval(std::make_pair(t0, t1)), _n0(n0), _n1(n1) {}

  int t0() const {return _interval.first;} // timestep
  int t1() const {return _interval.second;}
  int n0() const {return _n0;}
  int n1() const {return _n1;}

  ftkInterval GetInterval() const {return _interval;}
  void SetInterval(const ftkInterval &i) {_interval = i;}
 
public: // sparse matrix access
  inline void set(int i, int j, int val=1); 
  inline int get(int i, int j) const;

public: // events
  inline void detectEvents();
  const std::vector<ftkEvent>& events() const {return _events;}

private:
  ftkInterval _interval;
  int _n0 = 0, _n1 = 0;

  // matrix
  std::map<int, std::map<int, int> > _matrix; // row, colum, value

  // events
  std::vector<ftkEvent> _events;
};


////// 
void ftkTransitionMatrix::detectEvents()
{
  const int n = n0() + n1(); // number of graph nodes

  std::set<int> unvisited; 
  for (int i=0; i<n; i++) 
    unvisited.insert(i);

  while (!unvisited.empty()) {
    ftkEvent e;
    std::vector<int> Q;
    Q.push_back(*unvisited.begin());

    // find the left and right hand sides of the graph nodes
    while (!Q.empty()) {
      int v = Q.back();
      Q.pop_back();
      unvisited.erase(v);
      if (v<n0()) e.lhs.insert(v);
      else e.rhs.insert(v-n0());

      if (v<n0()) { // left hand side
        for (int j=0; j<n1(); j++) 
          if (get(v, j)>0 && unvisited.find(j+n0()) != unvisited.end())
            Q.push_back(j+n0());
      } else { // right hand side
        for (int i=0; i<n0(); i++) 
          if (get(i, v-n0())>0 && unvisited.find(i) != unvisited.end())
            Q.push_back(i);
      }
    }

    if (e.lhs.size() == 1 && e.rhs.size() == 1) { // no event
      e.type = FTK_EVENT_NONE;
    } else if (e.lhs.size() == 0 && e.rhs.size() == 1) {
      e.type = FTK_EVENT_BIRTH;
    } else if (e.lhs.size() == 1 && e.rhs.size() == 0) {
      e.type = FTK_EVENT_DEATH;
    } else if (e.lhs.size() == 1 && e.rhs.size() == 2) {
      e.type = FTK_EVENT_SPLIT;
    } else if (e.lhs.size() == 2 && e.rhs.size() == 1) { 
      e.type = FTK_EVENT_MERGE;
    } else if (e.lhs.size() == 2 && e.rhs.size() == 2) { 
      e.type = FTK_EVENT_RECOMBINATION;
    } else {
      e.type = FTK_EVENT_COMPOUND;
    }

    _events.push_back(e);
  }
}

void ftkTransitionMatrix::set(int i, int j, int val)
{
  if (val == 0) {
    if (_matrix.find(i) != _matrix.end()) {
      std::map<int, int>& vec = _matrix.at(i);
      if (vec.find(j) != vec.end()) 
        vec.erase(j);
      if (vec.empty())
        _matrix.erase(i);
    }
  } else 
    _matrix[i][j] = val;
}
  
inline int ftkTransitionMatrix::get(int i, int j) const 
{
  if (_matrix.find(i) == _matrix.end()) return 0;
  else {
    const std::map<int, int>& vec = _matrix.at(i);
    if (vec.find(j) == vec.cend()) return 0;
    else return vec.at(j);
  }
}

#endif
