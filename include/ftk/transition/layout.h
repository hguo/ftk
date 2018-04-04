#ifndef _FTK_LAYOUT_H
#define _FTK_LAYOUT_H

#include "ftk/transition/transition.h"
#include <fstream>

namespace ftk {

struct TransitionLayout {
  static void generateDotFile(const Transition&, const std::string& filename);
};

void TransitionLayout::generateDotFile(const Transition& tr, const std::string& filename)
{
  using namespace std;
  ofstream ofs(filename.c_str());
  if (!ofs.is_open()) return;

  auto node2str = [](int i, int j) {
    std::stringstream ss; 
    ss << i << "." << j;
    return ss.str();
  };

  ofs << "digraph {" << endl;
  ofs << "ratio = compress;" << endl;
  ofs << "rankdir = LR;" << endl;
  ofs << "ranksep =\"1.0 equally\";" << endl;
  ofs << "node [shape=circle];" << endl;
  // ofs << "node [shape=point];" << endl;
#if 0 // FIXME TODO
#if 1

  for (const auto &kv : _matrices) {
    const TransitionMatrix &tm = kv.second;
    for (int i=0; i<tm.n0(); i++) {
      for (int j=0; j<tm.n1(); j++) {
        int weight = 1;
        if (tm.rowsum(i) == 1 && tm.colsum(j) == 1) weight = 1000;

        if (tm.get(i, j)) {
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

  for (int t=_ts; t<_ts+_tl-1; t++) {
    const VortexTransitionMatrix &tm = Matrix(t); 
  }
#else
  // iteration over sequences
  for (int i=0; i<_seqs.size(); i++) {
    const VortexSequence &seq = _seqs[i];
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


}

#endif
