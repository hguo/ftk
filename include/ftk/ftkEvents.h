#ifndef _FEATURE_EVENTS_H
#define _FEATURE_EVENTS_H

#include <vector>
#include "common/Interval.h"

enum {
  FEATURE_EVENT_DUMMY = 0,
  FEATURE_EVENT_BIRTH = 1,
  FEATURE_EVENT_DEATH = 2,
  FEATURE_EVENT_MERGE = 3,
  FEATURE_EVENT_SPLIT = 4,
  FEATURE_EVENT_RECOMBINATION = 5, 
  FEATURE_EVENT_COMPOUND = 6
};

struct FeatureEvent {
  int if0, if1;
  int type;
  std::set<int> lhs, rhs; // local ids.
  // std::vector<int> lhs_gids, rhs_gids;

  std::map<int, float> dist; // key=timestep, val=dist

  static const char* TypeToString(int e) {
    static const char* strs[7] = {
      "dummy", "birth", "death", "merge", "split", "recombination", "compound"};
    return strs[e];
  }
};

#if 0
namespace diy {
  template <> struct Serialization<FeatureEvent> {
    static void save(diy::BinaryBuffer& bb, const FeatureEvent& m) {
      diy::save(bb, m.if0);
      diy::save(bb, m.if1);
      diy::save(bb, m.type);
      diy::save(bb, m.lhs);
      diy::save(bb, m.rhs);
    }

    static void load(diy::BinaryBuffer&bb, FeatureEvent& m) {
      diy::load(bb, m.if0);
      diy::load(bb, m.if1);
      diy::load(bb, m.type);
      diy::load(bb, m.lhs);
      diy::load(bb, m.rhs);
    }
  };
}
#endif

#endif
