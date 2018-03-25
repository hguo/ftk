#ifndef _FEATURE_TRANSITION_H
#define _FEATURE_TRANSITION_H

#include "common/FeatureTransitionMatrix.h"
#include "common/FeatureSequence.h"
#include <utility>
#include <mutex>

#if WITH_ROCKSDB
#include <rocksdb/db.h>
#endif

class FeatureTransition 
{
  // friend class diy::Serialization<FeatureTransition>;
public:
  FeatureTransition();
  ~FeatureTransition();
  
  // int ts() const {return _ts;}
  // int tl() const {return _tl;}

#ifdef WITH_ROCKSDB
  bool LoadFromDB(rocksdb::DB*);
#endif

  void LoadFromFile(const std::string &dataname, int ts, int tl);
  void SaveToDotFile(const std::string &filename) const;

  FeatureTransitionMatrix& Matrix(Interval intervals);
  void AddMatrix(const FeatureTransitionMatrix& m);
  int Transition(int t, int i, int j) const;
  // const std::map<int, FeatureTransitionMatrix>& Matrices() const {return _matrices;}
  std::map<Interval, FeatureTransitionMatrix>& Matrices() {return _matrices;}

  void ConstructSequence();
  void PrintSequence() const;
  void SequenceGraphColoring();
  void SequenceColor(int gid, unsigned char &r, unsigned char &g, unsigned char &b) const;

  int lvid2gvid(int frame, int lvid) const;
  int gvid2lvid(int frame, int gvid) const;

  int MaxNVorticesPerFrame() const {return _max_nvortices_per_frame;}
  int NVortices(int frame) const;

  const std::vector<struct FeatureSequence> Sequences() const {return _seqs;}
  void RandomColorSchemes();
  
  const std::vector<struct FeatureEvent>& Events() const {return _events;}

  int TimestepToFrame(int timestep) const {return _frames[timestep];} // confusing.  TODO: change func name
  int Frame(int i) const {return _frames[i];}
  int NTimesteps() const {return _frames.size();}
  void SetFrames(const std::vector<int> frames) {_frames = frames;}
  const std::vector<int>& Frames() const {return _frames;}

private:
  int NewFeatureSequence(int its);
  std::string NodeToString(int i, int j) const;

private:
  // int _ts, _tl;
  std::vector<int> _frames; // frame IDs
  std::map<Interval, FeatureTransitionMatrix> _matrices;
  std::vector<struct FeatureSequence> _seqs;
  std::map<std::pair<int, int>, int> _seqmap; // <time, lid>, gid
  std::map<std::pair<int, int>, int> _invseqmap; // <time, gid>, lid
  std::map<int, int> _nvortices_per_frame;
  int _max_nvortices_per_frame;

  std::vector<struct FeatureEvent> _events;

  std::mutex _mutex;
};

#if 0
/////////
namespace diy {
  template <> struct Serialization<FeatureTransition> {
    static void save(diy::BinaryBuffer& bb, const FeatureTransition& m) {
      diy::save(bb, m._frames);
      diy::save(bb, m._matrices);
      diy::save(bb, m._seqs);
      diy::save(bb, m._seqmap);
      diy::save(bb, m._invseqmap);
      diy::save(bb, m._nvortices_per_frame);
      diy::save(bb, m._max_nvortices_per_frame);
      diy::save(bb, m._events);
    }

    static void load(diy::BinaryBuffer&bb, FeatureTransition& m) {
      diy::load(bb, m._frames);
      diy::load(bb, m._matrices);
      diy::load(bb, m._seqs);
      diy::load(bb, m._seqmap);
      diy::load(bb, m._invseqmap);
      diy::load(bb, m._nvortices_per_frame);
      diy::load(bb, m._max_nvortices_per_frame);
      diy::load(bb, m._events);
    }
  };
}
#endif

#endif
