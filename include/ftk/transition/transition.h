#ifndef _FTK_TRANSITION_H
#define _FTK_TRANSITION_H

#include "ftk/transition/transitionMatrix.h"
#include "ftk/transition/sequence.h"
#include "ftk/storage/storage.h"
#include <utility>
#include <mutex>

class ftkTransition 
{
public:
  ftkTransition() {};
  ~ftkTransition() {};
  
  // int ts() const {return _ts;}
  // int tl() const {return _tl;}

  void LoadFromFile(const std::string &dataname, int ts, int tl) {} // FIXME
  // void SaveToDotFile(const std::string &filename) const;

  ftkTransitionMatrix& Matrix(const ftkInterval &I) {return _matrices[I];}
  inline void AddMatrix(const ftkTransitionMatrix& m);
  inline int Transition(int t, int i, int j) const;
  // const std::map<int, ftkTransitionMatrix>& Matrices() const {return _matrices;}
  std::map<ftkInterval, ftkTransitionMatrix>& Matrices() {return _matrices;}

  inline void ConstructSequence();
  inline void PrintSequence() const;
  // void SequenceGraphColoring();
  inline void SequenceColor(int gid, unsigned char &r, unsigned char &g, unsigned char &b) const;

  inline int lvid2gvid(int frame, int lvid) const;
  inline int gvid2lvid(int frame, int gvid) const;

  int MaxNVorticesPerFrame() const {return _max_nvortices_per_frame;}
  int NVortices(int frame) const;

  const std::vector<struct ftkSequence> Sequences() const {return _seqs;}
  void RandomColorSchemes();
  
  const std::vector<struct ftkEvent>& Events() const {return _events;}

  int TimestepToFrame(int timestep) const {return _frames[timestep];} // confusing.  TODO: change func name
  int Frame(int i) const {return _frames[i];}
  int NTimesteps() const {return _frames.size();}
  void SetFrames(const std::vector<int> frames) {_frames = frames;}
  const std::vector<int>& Frames() const {return _frames;}

private:
  inline int NewSequence(int its);
  std::string NodeToString(int i, int j) const;

private:
  // int _ts, _tl;
  std::vector<int> _frames; // frame IDs
  std::map<ftkInterval, ftkTransitionMatrix> _matrices;
  std::vector<struct ftkSequence> _seqs;
  std::map<std::pair<int, int>, int> _seqmap; // <time, lid>, gid
  std::map<std::pair<int, int>, int> _invseqmap; // <time, gid>, lid
  std::map<int, int> _nvortices_per_frame;
  int _max_nvortices_per_frame;

  std::vector<struct ftkEvent> _events;

  std::mutex _mutex;
};


//// impl
void ftkTransition::ConstructSequence()
{
  for (int i=0; i<_frames.size()-1; i++) {
    ftkInterval I(_frames[i], _frames[i+1]);
    // fprintf(stderr, "processing interval {%d, %d}\n", I.first, I.second);

    ftkTransitionMatrix &tm = Matrix(I);
    assert(tm.Valid());

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

void ftkTransition::AddMatrix(const ftkTransitionMatrix& m)
{
  if (!m.Valid()) return;
  
  std::unique_lock<std::mutex> lock(_mutex);
  _matrices[m.GetInterval()] = m;
}

void ftkTransition::SequenceColor(int gid, unsigned char &r, unsigned char &g, unsigned char &b) const
{
  r = _seqs[gid].r;
  g = _seqs[gid].g;
  b = _seqs[gid].b;
}

#endif
