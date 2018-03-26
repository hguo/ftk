#ifndef _FTK_TRANSITION_H
#define _FTK_TRANSITION_H

#include "common/ftkTransitionMatrix.h"
#include "common/ftkSequence.h"
#include <utility>
#include <mutex>

#if WITH_ROCKSDB
#include <rocksdb/db.h>
#endif

class ftkTransition 
{
  // friend class diy::Serialization<ftkTransition>;
public:
  ftkTransition();
  ~ftkTransition();
  
  // int ts() const {return _ts;}
  // int tl() const {return _tl;}

#ifdef WITH_ROCKSDB
  bool LoadFromDB(rocksdb::DB*);
#endif

  void LoadFromFile(const std::string &dataname, int ts, int tl);
  void SaveToDotFile(const std::string &filename) const;

  ftkTransitionMatrix& Matrix(Interval intervals);
  void AddMatrix(const ftkTransitionMatrix& m);
  int Transition(int t, int i, int j) const;
  // const std::map<int, ftkTransitionMatrix>& Matrices() const {return _matrices;}
  std::map<Interval, ftkTransitionMatrix>& Matrices() {return _matrices;}

  void ConstructSequence();
  void PrintSequence() const;
  void SequenceGraphColoring();
  void SequenceColor(int gid, unsigned char &r, unsigned char &g, unsigned char &b) const;

  int lvid2gvid(int frame, int lvid) const;
  int gvid2lvid(int frame, int gvid) const;

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
  int NewftkSequence(int its);
  std::string NodeToString(int i, int j) const;

private:
  // int _ts, _tl;
  std::vector<int> _frames; // frame IDs
  std::map<Interval, ftkTransitionMatrix> _matrices;
  std::vector<struct ftkSequence> _seqs;
  std::map<std::pair<int, int>, int> _seqmap; // <time, lid>, gid
  std::map<std::pair<int, int>, int> _invseqmap; // <time, gid>, lid
  std::map<int, int> _nvortices_per_frame;
  int _max_nvortices_per_frame;

  std::vector<struct ftkEvent> _events;

  std::mutex _mutex;
};

#endif
