#ifndef _FTK_EVENTS_H
#define _FTK_EVENTS_H

#include <vector>
#include <set>
#include "ftk/base/object.h"
#include "ftk/transition/interval.h"
#include "ftk/external/json.hh"

enum {
  FTK_EVENT_NONE = 0, // no event
  FTK_EVENT_BIRTH = 1,
  FTK_EVENT_DEATH = 2,
  FTK_EVENT_MERGE = 3,
  FTK_EVENT_SPLIT = 4,
  FTK_EVENT_RECOMBINATION = 5, 
  FTK_EVENT_COMPOUND = 6
};

struct ftkEvent {
  // int if0, if1;
  int type;
  std::set<int> lhs, rhs; // local ids.
  // std::vector<int> lhs_gids, rhs_gids;

  static const char* TypeToString(int e) {
    static const char* strs[7] = {
      "none", "birth", "death", "merge", "split", "recombination", "compound"};
    return strs[e];
  }
};

#endif
