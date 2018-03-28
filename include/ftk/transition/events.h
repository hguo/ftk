#ifndef _FTK_EVENTS_H
#define _FTK_EVENTS_H

#include <vector>
#include <set>
#include "ftk/base/object.h"
#include "ftk/transition/interval.h"
#include "ftk/external/json.hh"

enum {
  FTK_EVENT_DUMMY = 0,
  FTK_EVENT_BIRTH = 1,
  FTK_EVENT_DEATH = 2,
  FTK_EVENT_MERGE = 3,
  FTK_EVENT_SPLIT = 4,
  FTK_EVENT_RECOMBINATION = 5, 
  FTK_EVENT_COMPOUND = 6
};

struct ftkEvent : public ftkObject {
  int if0, if1;
  int type;
  std::set<int> lhs, rhs; // local ids.
  // std::vector<int> lhs_gids, rhs_gids;

  static const char* TypeToString(int e) {
    static const char* strs[7] = {
      "dummy", "birth", "death", "merge", "split", "recombination", "compound"};
    return strs[e];
  }


  using json = nlohmann::json;

  json toJson() const {
    json j;
    
    j["if0"] = if0; 
    j["if1"] = if1;
    j["type"] = type;
    j["lhs"] = lhs;
    j["rhs"] = rhs;

    return j;
  }

  void fromJson(json j) {
    if0 = j["if0"];
    if1 = j["if1"];
    type = j["type"];
    lhs = j["lhs"].get<std::set<int> >();
    rhs = j["rhs"].get<std::set<int> >();
  }
};

#endif
