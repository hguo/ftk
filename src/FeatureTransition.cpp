#include "FeatureTransition.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <set>
#include <cassert>
#include <cstring>
// #include "common/diy-ext.hpp"
#include "RandomColor.h"
#include "GraphColor.h"

// mem alloc for 2D arrays
#define malloc2D(name, xDim, yDim, type) do {               \
    name = (type **)malloc(xDim * sizeof(type *));          \
    assert(name != NULL);                                   \
    name[0] = (type *)malloc(xDim * yDim * sizeof(type));   \
    assert(name[0] != NULL);                                \
    for (size_t i = 1; i < xDim; i++)                       \
        name[i] = name[i-1] + yDim;                         \
} while (0)

FeatureTransition::FeatureTransition()
{
}

FeatureTransition::~FeatureTransition()
{
}

#if WITH_ROCKSDB
bool FeatureTransition::LoadFromDB(rocksdb::DB* db)
{
  std::string buf;
  rocksdb::Status s; 

  s = db->Get(rocksdb::ReadOptions(), "trans", &buf);
  if (s.ok()) {
    diy::unserialize(buf, *this);
  } else {
    s = db->Get(rocksdb::ReadOptions(), "f", &buf);
    if (!s.ok()) return false;

    diy::unserialize(buf, _frames);
    const int nframes = _frames.size();

    fprintf(stderr, "nframes=%d\n", nframes);

    for (int i=0; i<nframes-1; i++) {
      std::stringstream ss;
      ss << "m." << _frames[i] << "." << _frames[i+1];
      rocksdb::Status s = db->Get(rocksdb::ReadOptions(), ss.str(), &buf);
      if (!s.ok()) fprintf(stderr, "Key not found, %s\n", ss.str().c_str());

      FeatureTransitionMatrix mat;
      diy::unserialize(buf, mat);
      AddMatrix(mat);
    }

    ConstructSequence();
    diy::serialize(*this, buf);
    db->Put(rocksdb::WriteOptions(), "trans", buf);
  }

  return true;
}
#endif

void FeatureTransition::LoadFromFile(const std::string& dataname, int ts, int tl)
{
#if 0 // FIXME
  _ts = ts;
  _tl = tl;

  for (int i=ts; i<ts+tl-1; i++) {
    std::stringstream ss;
    ss << dataname << ".match." << i << "." << i+1;

    FeatureTransitionMatrix tm;
    if (!tm.LoadFromFile(ss.str())) {
      fprintf(stderr, "cannot open file %s\n", ss.str().c_str());
      tm.SetToDummy();
    }
    
    _matrices.insert(std::make_pair(i, tm));
    _nvortices_per_frame[i] = tm.n0(); 
    _nvortices_per_frame[i+1] = tm.n1();
  }

  _max_nvortices_per_frame = 0;
  for (std::map<int, int>::iterator it = _nvortices_per_frame.begin(); it != _nvortices_per_frame.end(); it ++) {
    _max_nvortices_per_frame = std::max(_max_nvortices_per_frame, it->second);
  }
  // fprintf(stderr, "max_nvortices_per_frame=%d\n", _max_nvortices_per_frame);
#endif
}

#define HUMAN_READABLE 0
std::string FeatureTransition::NodeToString(int i, int j) const
{
  using namespace std;
  stringstream ss;
#if HUMAN_READABLE
  ss << i << "." << j;
#else
  ss << i*_max_nvortices_per_frame + j;
#endif
  return ss.str();
}

void FeatureTransition::SaveToDotFile(const std::string& filename) const
{
  using namespace std;
  ofstream ofs(filename.c_str());
  if (!ofs.is_open()) return;

  ofs << "digraph {" << endl;
  ofs << "ratio = compress;" << endl;
  ofs << "rankdir = LR;" << endl;
  ofs << "ranksep =\"1.0 equally\";" << endl;
  ofs << "node [shape=circle];" << endl;
  // ofs << "node [shape=point];" << endl;
#if 0 // FIXME TODO
#if 1
  for (int t=_ts; t<_ts+_tl-1; t++) {
    const FeatureTransitionMatrix &tm = Matrix(t); 
    for (int i=0; i<tm.n0(); i++) {
      for (int j=0; j<tm.n1(); j++) {
        int weight = 1;
        if (tm.rowsum(i) == 1 && tm.colsum(j) == 1) weight = 1000;

        if (tm(i, j)) {
          ofs << NodeToString(t, i) << "->" 
              << NodeToString(t+1, j)
              << " [weight = " << weight << "];" << endl;
        }
      }
    }
    
    ofs << "{ rank=same; ";
    for (int i=0; i<tm.n0(); i++) {
      if (i<tm.n0()-1) ofs << NodeToString(t, i) << ", ";
      else ofs << NodeToString(t, i) << "}" << endl;
    }
  }
#else
  // iteration over sequences
  for (int i=0; i<_seqs.size(); i++) {
    const FeatureSequence &seq = _seqs[i];
    for (int k=0; k<seq.lids.size(); k++) {
      const int t = seq.ts + k;
      const int weight = seq.tl;
      if (k<seq.lids.size()-1) 
        ofs << NodeToString(t, seq.lids[k]) << "->";
      else 
        ofs << NodeToString(t, seq.lids[k]) 
            << " [weight = " << weight << "];" << endl;
    }
  }
  // ranks
  for (int t=_ts; t<_ts+_tl-1; t++) {
    std::map<int, int>::const_iterator it = _nvortices_per_frame.find(t); 
    const int n = it->second;
    ofs << "{ rank=same; ";
    for (int i=0; i<n; i++) {
      if (i<n-1) ofs << NodeToString(t, i) << ", ";
      else ofs << NodeToString(t, i) << " }" << endl;
    }
  }
  // subgraphs
  for (int i=0; i<_events.size(); i++) {
    int t = _events[i].t;
    ofs << "subgraph {";
    for (std::set<int>::iterator it0 = _events[i].lhs.begin(); it0 != _events[i].lhs.end(); it0 ++) 
      for (std::set<int>::iterator it1 = _events[i].rhs.begin(); it1 != _events[i].rhs.end(); it1 ++) {
        ofs << NodeToString(t, *it0) << "->" 
            << NodeToString(t+1, *it1) << ";" << endl;
      }
    ofs << "};" << endl;
  }
#endif
#endif
  // node colors
#if 0 // TODO FIXME
  for (int t=_ts; t<_ts+_tl; t++) {
    std::map<int, int>::const_iterator it = _nvortices_per_frame.find(t); 
    const int n = it->second;
    for (int k=0; k<n; k++) {
      const int nc = 6;
      int vid = lvid2gvid(t, k);
      int c = vid % nc;
      std::string color;
      
      if (c == 0) color = "blue";
      else if (c == 1) color = "green";
      else if (c == 2) color = "cyan";
      else if (c == 3) color = "red";
      else if (c == 4) color = "purple";
      else if (c == 5) color = "yellow";

#if HUMAN_READABLE
      ofs << t << "." << k 
          << " [style=filled, fillcolor=" << color << "];" << endl;
#else
      ofs << t*_max_nvortices_per_frame+k
          << " [style=filled, fillcolor=" << color << "];" << endl;
#endif
    }
  }
  ofs << "}" << endl;
  ofs.close();
#endif
}

FeatureTransitionMatrix& FeatureTransition::Matrix(Interval I)
{
  return _matrices[I];
}

void FeatureTransition::AddMatrix(const FeatureTransitionMatrix& m)
{
  if (!m.Valid()) return;
  
  std::unique_lock<std::mutex> lock(_mutex);
  _matrices[m.GetInterval()] = m;
}

int FeatureTransition::NewFeatureSequence(int its)
{
  FeatureSequence vs;
  vs.its = its;
  vs.itl = 0;
  // vs.lhs_event = vs.rhs_event = FEATURE_EVENT_DUMMY;
  _seqs.push_back(vs);
  return _seqs.size() - 1;
}

int FeatureTransition::lvid2gvid(int t, int lid) const
{
  std::pair<int, int> key = std::make_pair(t, lid);
  std::map<std::pair<int,int>,int>::const_iterator it = _seqmap.find(key);
  if (it == _seqmap.end())
    return -1;
  else 
    return it->second;
}

int FeatureTransition::gvid2lvid(int frame, int gvid) const
{
  std::pair<int, int> key = std::make_pair(frame, gvid);
  std::map<std::pair<int, int>,int>::const_iterator it = _invseqmap.find(key);
  if (it == _invseqmap.end())
    return -1;
  else 
    return it->second;
}

void FeatureTransition::SequenceColor(int gid, unsigned char &r, unsigned char &g, unsigned char &b) const
{
  r = _seqs[gid].r;
  g = _seqs[gid].g;
  b = _seqs[gid].b;
}

void FeatureTransition::ConstructSequence()
{
  for (int i=0; i<_frames.size()-1; i++) {
    Interval I(_frames[i], _frames[i+1]);
    // fprintf(stderr, "processing interval {%d, %d}\n", I.first, I.second);

    FeatureTransitionMatrix &tm = Matrix(I);
    assert(tm.Valid());

    if (i == 0) { // initial
      for (int k=0; k<tm.n0(); k++) {
        int gid = NewFeatureSequence(i);
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
          int gid = NewFeatureSequence(i+1);
          _seqs[gid].itl ++;
          _seqs[gid].lids.push_back(r);
          _seqmap[std::make_pair(i+1, r)] = gid;
          _invseqmap[std::make_pair(i+1, gid)] = r;
        }
      }

      // build events
      // if (event >= FEATURE_EVENT_MERGE) {
      if (event > FEATURE_EVENT_DUMMY) {
        FeatureEvent e;
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
  SequenceGraphColoring(); 
}

void FeatureTransition::PrintSequence() const
{
  for (int i=0; i<_events.size(); i++) {
    const FeatureEvent& e = _events[i];
    std::stringstream ss;
    ss << "interval={" << _frames[e.if0] << ", " << _frames[e.if1] << "}, ";
    ss << "type=" << FeatureEvent::TypeToString(e.type) << ", ";
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


#if 0
      if (event >= FEATURE_EVENT_MERGE) { // FIXME: other events
        for (int u=0; u<lhs.size(); u++) {
          const int lgid = _seqmap[std::make_tuple(i, lhs[u])];
          _seqs[lgid].rhs_event = event;
          _seqs[lgid].rhs_gids = rhs;
        }
      }
#endif


#if 0
      int cnt=0;
      fprintf(stderr, "lhs={");
      for (std::set<int>::iterator it=lhs.begin(); it!=lhs.end(); it++) {
        if (cnt<lhs.size()-1) fprintf(stderr, "%d, ", *it);
        else fprintf(stderr, "%d}, ", *it);
        cnt ++;
      }
      cnt=0;
      fprintf(stderr, "rhs={");
      for (std::set<int>::iterator it=rhs.begin(); it!=rhs.end(); it++) {
        if (cnt<rhs.size()-1) fprintf(stderr, "%d, ", *it);
        else fprintf(stderr, "%d}\n", *it);
        cnt ++;
      }
#endif

void FeatureTransition::RandomColorSchemes()
{
  std::vector<unsigned char> colors;
  generate_random_colors(_seqs.size(), colors);

  for (int i=0; i<_seqs.size(); i++) {
    _seqs[i].r = colors[i*3];
    _seqs[i].g = colors[i*3+1];
    _seqs[i].b = colors[i*3+2];
  }
}

void FeatureTransition::SequenceGraphColoring()
{
  using namespace std;

  // 1. construct graph
  const int n = _seqs.size();
  bool **M;
  malloc2D(M, n, n, bool);
  memset(M[0], 0, sizeof(bool)*n*n);

  // 1.1 concurrent vortices
  map<int, std::set<int> > time_slots;
  for (int i=0; i<n; i++) {
    for (int j=0; j<_seqs[i].itl; j++) {
      time_slots[_seqs[i].its+j].insert(i);
    }
  }
  for (map<int, std::set<int> >::iterator it=time_slots.begin(); it!=time_slots.end(); it++) {
    const set<int> &concurrents = it->second;
    for (set<int>::iterator it0=concurrents.begin(); it0!=concurrents.end(); it0++)
      for (set<int>::iterator it1=concurrents.begin(); it1!=concurrents.end(); it1++)
        if (it0!=it1) {
          M[*it0][*it1] = M[*it1][*it0] = true;
        }
  }

  // 1.2 events
  for (int i=0; i<n; i++) {
    int t = _seqs[i].its + _seqs[i].itl - 1;
    if (t>=_frames.size()-1) continue; 
    int lhs_lid = _seqs[i].lids.back();
    Interval interval = std::make_pair(_frames[t], _frames[t+1]);
    FeatureTransitionMatrix &mat = _matrices[interval];
    for (int k=0; k<mat.n1(); k++) {
      if (mat(lhs_lid, k)) {
        int rhs_lid = k;
        int rgid = _seqmap[std::make_pair(t+1, k)];
        M[i][rgid] = M[rgid][i] = true;
      }
    }
  }

  // 2. graph coloring
  int *cids = (int*)malloc(sizeof(int)*n);
  int nc = welsh_powell(n, M, cids);

  // 3. generate colors
  // fprintf(stderr, "#color=%d\n", nc);
  vector<unsigned char> colors;
  generate_random_colors(nc, colors);

  for (int i=0; i<_seqs.size(); i++) {
    int c = cids[i];
    // fprintf(stderr, "seq=%d, c=%d\n", i, c);
    _seqs[i].r = colors[c*3];
    _seqs[i].g = colors[c*3+1];
    _seqs[i].b = colors[c*3+2];
  }

  // 4. release memory
  free(M[0]);
  free(M);
  free(cids);
}

int FeatureTransition::NVortices(int frame) const
{
  std::map<int, int>::const_iterator it = _nvortices_per_frame.find(frame);

  if (it != _nvortices_per_frame.end())
    return it->second;
  else 
    return 0;
}
