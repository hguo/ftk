#ifndef _FTK_TRANSITION_H
#define _FTK_TRANSITION_H

#include "ftk/transition/transitionMatrix.h"
#include "ftk/transition/event.h"
#include "ftk/storage/storage.h"
#include <utility>
#include <mutex>

class ftkTransition 
{
public:
  ftkTransition() {};
  ~ftkTransition() {};

  inline void addInterval(int t0, int n0, int t1, int n1);
  inline void addTransition(int t0, int lid0, int t1, int lid1);

  inline void relabel();

  inline int lid2gid(int t, int lid) const;
  inline int gid2lid(int t, int gid) const;

  inline void printComponents() const;

private:
  std::map<int, int> _intervals; // <time, next time>
  std::map<int, ftkTransitionMatrix> _matrices; // <time, matrix>
  std::vector<ftkEvent> _events;
  
  std::map<int, std::map<int, int> > _components; // gid, <t, lid>
  std::map<std::pair<int, int>, int> _labels; // <t, lid>, gid

  mutable std::mutex _mutex;
};


//// impl
void ftkTransition::addInterval(int t0, int n0, int t1, int n1)
{
  std::unique_lock<std::mutex> lock(_mutex);
  _intervals[t0] = t1;
  _matrices[t0] = ftkTransitionMatrix(n0, n1);
}

void ftkTransition::addTransition(int t0, int lid0, int t1, int lid1)
{
  std::unique_lock<std::mutex> lock(_mutex);
  _matrices[t0].set(lid0, lid1);
}
  
int ftkTransition::lid2gid(int t, int lid) const 
{
  std::unique_lock<std::mutex> lock(_mutex);
  auto key = std::make_pair(t, lid);
  if (_labels.find(key) == _labels.cend()) return -1; // not found
  else return _labels.at(key);
}

int ftkTransition::gid2lid(int t, int gid) const 
{
  std::unique_lock<std::mutex> lock(_mutex);
  if (_components.find(gid) == _components.cend()) return -1;
  else {
    const auto &component = _components.at(gid);
    if (component.find(t) == component.cend()) return -1; 
    else return component.at(t);
  }
}

void ftkTransition::relabel()
{
  std::unique_lock<std::mutex> lock(_mutex);
  
  bool first = true;
  int nComponents = 0;

  _components.clear();
  _labels.clear();

  for (auto &kv : _matrices) {
    const int t = kv.first;
    const int t1 = _intervals[t];
    ftkTransitionMatrix& m = kv.second;
  
    if (first) {
      first = false;
      for (int k=0; k<m.n0(); k++) {
        int gid = nComponents ++; // new component
        _components[gid][t] = k;
        _labels[std::make_pair(t, k)] = gid;
      }
    }

    if (m.events().size() == 0) m.detectEvents();
    for (const auto &e : m.events()) {
      if (e.lhs.size() == 1 && e.rhs.size() == 1) { // no critical event
        int l = *e.lhs.begin(), r = *e.rhs.begin();
        int gid = _labels[std::make_pair(t, l)];
        _components[gid][t1] = r;
        _labels[std::make_pair(t1, r)] = gid;
      } else { // event
        for (const auto &r : e.rhs) {
          int gid = nComponents ++; // new component
          _components[gid][t1] = r;
          _labels[std::make_pair(t1, r)] = gid;
        }
      }
    }
  }

  // printComponents();
}

void ftkTransition::printComponents() const 
{
  std::unique_lock<std::mutex> lock(_mutex);
  
  for (const auto &kv : _components) {
    fprintf(stderr, "componentId=%d\n", kv.first); 
    for (const auto &kv1 : kv.second) {
      fprintf(stderr, "t=%d, lid=%d\n", kv1.first, kv1.second);
    }
  }
}

#endif
