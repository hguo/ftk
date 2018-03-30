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

  std::mutex _mutex;
};


//// impl
void ftkTransition::addInterval(int t0, int n0, int t1, int n1)
{
  _intervals[t0] = t1;
  _matrices[t0] = ftkTransitionMatrix(t0, n0, t1, n1);
}

void ftkTransition::addTransition(int t0, int lid0, int t1, int lid1)
{
  _matrices[t0].set(lid0, lid1);
}
  
int ftkTransition::lid2gid(int t, int lid) const 
{
  auto key = std::make_pair(t, lid);
  if (_labels.find(key) == _labels.cend()) return -1; // not found
  else return _labels.at(key);
}

int ftkTransition::gid2lid(int t, int gid) const 
{
  if (_components.find(gid) == _components.cend()) return -1;
  else {
    const auto &component = _components.at(gid);
    if (component.find(t) == component.cend()) return -1; 
    else return component.at(t);
  }
}

void ftkTransition::relabel()
{
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
  for (const auto &kv : _components) {
    fprintf(stderr, "componentId=%d\n", kv.first); 
    for (const auto &kv1 : kv.second) {
      fprintf(stderr, "t=%d, lid=%d\n", kv1.first, kv1.second);
    }
  }
}

#if 0
void ftkTransition::ConstructSequence()
{
  for (int i=0; i<_frames.size()-1; i++) {
    ftkInterval I(_frames[i], _frames[i+1]);
    ftkTransitionMatrix &tm = Matrix(I);

    if (i == 0) { // initial
      for (int k=0; k<tm.n0(); k++) {
        int gid = NewSequence(i);
        _seqs[gid].itl ++;
        _seqs[gid].lids.push_back(k);
        _seqmap[std::make_pair(i, k)] = gid;
        _invseqmap[std::make_pair(i, gid)] = k;
      }
    }

    for (int k=0; k<tm.NModules(); k++) {
      int event;
      std::set<int> lhs, rhs;
      tm.GetModule(k, lhs, rhs, event);

      if (lhs.size() == 1 && rhs.size() == 1) { // ordinary case
        int l = *lhs.begin(), r = *rhs.begin();
        int gid = _seqmap[std::make_pair(i, l)];
        _seqs[gid].itl ++;
        _seqs[gid].lids.push_back(r);
        _seqmap[std::make_pair(i+1, r)] = gid;
        _invseqmap[std::make_pair(i+1, gid)] = r;
      } else { // some events, need re-ID
        for (std::set<int>::iterator it=rhs.begin(); it!=rhs.end(); it++) {
          int r = *it; 
          int gid = NewSequence(i+1);
          _seqs[gid].itl ++;
          _seqs[gid].lids.push_back(r);
          _seqmap[std::make_pair(i+1, r)] = gid;
          _invseqmap[std::make_pair(i+1, gid)] = r;
        }
      }

      // build events
      if (event > FTK_EVENT_DUMMY) {
        ftkEvent e;
        e.if0 = i; 
        e.if1 = i+1;
        e.type = event;
        e.lhs = lhs;
        e.rhs = rhs;
        _events.push_back(e);
      }
    }
  }

  // RandomColorSchemes();
  // SequenceGraphColoring(); 
}

int ftkTransition::NewSequence(int its)
{
  ftkSequence vs;
  vs.its = its;
  vs.itl = 0;
  // vs.lhs_event = vs.rhs_event = VORTEX_EVENT_DUMMY;
  _seqs.push_back(vs);
  return _seqs.size() - 1;
}
 
int ftkTransition::lvid2gvid(int t, int lid) const
{
  std::pair<int, int> key = std::make_pair(t, lid);
  std::map<std::pair<int,int>,int>::const_iterator it = _seqmap.find(key);
  if (it == _seqmap.end())
    return -1;
  else 
    return it->second;
}

int ftkTransition::gvid2lvid(int frame, int gvid) const
{
  std::pair<int, int> key = std::make_pair(frame, gvid);
  std::map<std::pair<int, int>,int>::const_iterator it = _invseqmap.find(key);
  if (it == _invseqmap.end())
    return -1;
  else 
    return it->second;
}

void ftkTransition::PrintSequence() const
{
  for (int i=0; i<_events.size(); i++) {
    const ftkEvent& e = _events[i];
    std::stringstream ss;
    ss << "interval={" << _frames[e.if0] << ", " << _frames[e.if1] << "}, ";
    ss << "type=" << ftkEvent::TypeToString(e.type) << ", ";
    ss << "lhs={";

    int j = 0;
    if (e.lhs.empty()) ss << "}, "; 
    else 
      for (std::set<int>::iterator it = e.lhs.begin(); it != e.lhs.end(); it++, j++) {
        const int gvid = lvid2gvid(e.if0, *it);
        if (j<e.lhs.size()-1) 
          ss << gvid << ", ";
        else 
          ss << gvid << "}, ";
      }
    
    ss << "rhs={";
   
    j = 0;
    if (e.rhs.empty()) ss << "}";
    else 
      for (std::set<int>::iterator it = e.rhs.begin(); it != e.rhs.end(); it++, j++) {
        const int gvid = lvid2gvid(e.if1, *it);
        if (j<e.rhs.size()-1) 
          ss << gvid << ", ";
        else 
          ss << gvid << "}";
      }
    
    std::cout << ss.str() << std::endl;
  }
}

void ftkTransition::SequenceColor(int gid, unsigned char &r, unsigned char &g, unsigned char &b) const
{
  r = _seqs[gid].r;
  g = _seqs[gid].g;
  b = _seqs[gid].b;
}
#endif


#endif
