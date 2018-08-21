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
  friend class nlohmann::adl_serializer<ftk::TransitionMatrix>;

public:
  TransitionMatrix() {};
  TransitionMatrix(int t0, int n0, int t1, int n1) : _interval(std::make_pair(t0, t1)), _n0(n0), _n1(n1) {}

  std::pair<int, int> interval() const {return _interval;}
  int n0() const {return _n0;}
  int n1() const {return _n1;}
  int t0() const {return _interval.first;}
  int t1() const {return _interval.second;}

public: // sparse matrix access
  inline void set(int i, int j, int val=1); 
  inline int get(int i, int j) const;

  inline int colsum(int j) const;
  inline int rowsum(int i) const;

public: // events
  inline void detectEvents();
  const std::vector<Event<int, int> >& events() const {return _events;}

private:
  int _n0 = 0, _n1 = 0;
  std::pair<int, int> _interval;
  std::map<int, std::map<int, int> > _matrix; // row, colum, value

  // events
  std::vector<Event<int, int> > _events;
};


////// 
int TransitionMatrix::rowsum(int i) const
{
  if (_matrix.find(i) != _matrix.end()) {
    const std::map<int, int>& vec = _matrix.at(i);
    int sum = 0;
    for (const auto &kv : vec) {
      sum += kv.second;
    }
    return sum;
  } else return 0;
}

int TransitionMatrix::colsum(int j) const 
{
  int sum = 0;
  for (const auto &kv : _matrix) {
    if (kv.second.find(j) != kv.second.end()) 
      sum += kv.second.at(j);
  }
  return sum;
}

void TransitionMatrix::detectEvents()
{
  const int n = n0() + n1(); // number of graph nodes

  std::set<int> unvisited; 
  for (int i=0; i<n; i++) 
    unvisited.insert(i);

  while (!unvisited.empty()) {
    Event<int, int> e;
    e.interval = _interval;
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

// json
namespace nlohmann {
  template <>
  struct adl_serializer<ftk::TransitionMatrix> {
    static void to_json(json& j, const ftk::TransitionMatrix &m) {
      j["t0"] = m._interval.first;
      j["t1"] = m._interval.second;
      j["n0"] = m._n0;
      j["n1"] = m._n1;

      std::vector<int> d;
      for (const auto &kv : m._matrix) {
        const int i = kv.first;
        for (const auto &kv1 : kv.second) {
          const int j = kv1.first;
          // const int val = kv1.second;
          d.push_back(i);
          d.push_back(j);
          // d.push_back(val);
        }
      }
      j["transition"] = d;
    }

    static void from_json(const json& j, ftk::TransitionMatrix &m) {
      m._interval = std::make_pair<int, int>(j["t0"], j["t1"]);
      m._n0 = j["n0"];
      m._n1 = j["n1"];

      std::vector<int> d = j["transition"];
      for (int i=0; i<d.size()/2; i++)
        m.set(d[i*2], d[i*2+1]);
    }
  };
}

#endif
