#ifndef _FTK_LAYOUT_H
#define _FTK_LAYOUT_H

#include "ftk/transition/transition.h"
#include <fstream>

namespace ftk {

struct TransitionLayout {
  static void generateDotInputFile(const Transition&, const std::string& filename);
  static void parseDotOutputFile(const Transition&, const std::string& filename);
};

void TransitionLayout::generateDotInputFile(const Transition& tr, const std::string& filename)
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
  // ofs << "node [shape=circle];" << endl;
  ofs << "node [shape=point,width=0,height=0];" << endl;
  
  // nodes and links
  for (const auto &kv : tr.matrices()) {
    const TransitionMatrix &tm = kv.second;
    const Interval interval = kv.second.interval();

    for (int i=0; i<tm.n0(); i++) {
      for (int j=0; j<tm.n1(); j++) {
        int weight = 1;
        if (tm.rowsum(i) == 1 && tm.colsum(j) == 1) weight = 1000;

        if (tm.get(i, j)) {
          ofs << node2str(interval.first, i) << "->" 
              << node2str(interval.second, j)
              << " [weight = " << weight << "];" << endl;
        }
      }
    }
    
    ofs << "{ rank=same; ";
    for (int i=0; i<tm.n0(); i++) {
      if (i<tm.n0()-1) ofs << node2str(interval.first, i) << ", ";
      else ofs << node2str(interval.first, i) << "}" << endl;
    }
  }

  // node colors
  for (const auto& kv : tr.labels()) {
    int c = kv.second % 6;
    std::string color;
      
    if (c == 0) color = "blue";
    else if (c == 1) color = "green";
    else if (c == 2) color = "cyan";
    else if (c == 3) color = "red";
    else if (c == 4) color = "purple";
    else if (c == 5) color = "yellow";

    ofs << node2str(kv.first.first, kv.first.second) 
        << " [style=filled, fillcolor=" << color << "];" << endl;
  }
  
  ofs << "}" << endl;
  ofs.close();
}

void TransitionLayout::parseDotOutputFile(const Transition&, const std::string& filename)
{
  using namespace std;
  ifstream ifs(filename);
  if (!ifs.is_open()) return;

  float layout_width, layout_height;
  string str, label, type, shape, color, fillcolor;
  std::string id;
  float x, y, w, h;

  ifs >> str >> id /* dummy */ >> layout_width >> layout_height;
  fprintf(stderr, "layout_width=%f, layout_height=%f\n", layout_width, layout_height);
  while (1) {
    ifs >> str;
    if (str == "node") {
      ifs >> id >> x >> y >> w >> h >> label >> type >> shape >> color >> fillcolor;
      fprintf(stderr, "id=%s, x=%f, y=%f, w=%f, h=%f\n", id.c_str(), x, y, w, h);
    } else 
      break;
  }

  ifs.close();
}

}

#endif
