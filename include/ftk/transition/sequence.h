#ifndef _FTK_SEQUENCE_H
#define _FTK_SEQUENCE_H

#include <vector>
#include <map>
#include "ftk/transition/transition.h"

struct ftkSequence {
  // int ts, tl; // start and duration
  int its, itl;  // start and duration (index of frames)
  std::vector<int> lids; // local ids

  unsigned char r, g, b;

  // std::vector<int> lhs_gids, rhs_gids;
  // int lhs_event, rhs_event;

  std::string toJson() const {
    nlohmann::json j;
    
    j["its"] = its; 
    j["itl"] = itl;
    j["lids"] = lids;
    j["r"] = r; 
    j["g"] = g;
    j["b"] = b;

    return j.dump();
  }

  void fromJson(const std::string &str) {
    auto j = nlohmann::json::parse(str);
    its = j["its"];
    itl = j["itl"];
    r = j["r"];
    g = j["g"];
    b = j["b"];
  }
};

#endif
