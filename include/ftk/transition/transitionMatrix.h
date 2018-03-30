#ifndef _FTK_TRANSITION_MATRIX_H
#define _FTK_TRANSITION_MATRIX_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include "ftk/transition/interval.h"
#include "ftk/transition/event.h"

namespace ftk {

class TransitionMatrix {
public:
  TransitionMatrix() {};
  TransitionMatrix(int n0, int n1) : _n0(n0), _n1(n1) {}

  int n0() const {return _n0;}
  int n1() const {return _n1;}

public: // sparse matrix access
  inline void set(int i, int j, int val=1); 
  inline int get(int i, int j) const;

public: // events
  inline void detectEvents();
  const std::vector<Event>& events() const {return _events;}

private:
  int _n0 = 0, _n1 = 0;
  std::map<int, std::map<int, int> > _matrix; // row, colum, value

  // events
  std::vector<Event> _events;
};


////// 
void TransitionMatrix::detectEvents()
{
  const int n = n0() + n1(); // number of graph nodes

  std::set<int> unvisited; 
  for (int i=0; i<n; i++) 
    unvisited.insert(i);

  while (!unvisited.empty()) {
    Event e;
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

    _events.push_back(e);
  }
}

void TransitionMatrix::set(int i, int j, int val)
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
  
inline int TransitionMatrix::get(int i, int j) const 
{
  if (_matrix.find(i) == _matrix.end()) return 0;
  else {
    const std::map<int, int>& vec = _matrix.at(i);
    if (vec.find(j) == vec.cend()) return 0;
    else return vec.at(j);
  }
}

}

#endif
