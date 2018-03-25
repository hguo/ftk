#ifndef _FEATURE_SEQUENCE_H
#define _FEATURE_SEQUENCE_H

#include <vector>
#include <map>
#include "common/FeatureTransition.h"

struct FeatureSequence {
  // int ts, tl; // start and duration
  int its, itl;  // start and duration (index of frames)
  std::vector<int> lids; // local ids

  unsigned char r, g, b;

  // std::vector<int> lhs_gids, rhs_gids;
  // int lhs_event, rhs_event;
};

#if 0
namespace diy {
  template <> struct Serialization<FeatureSequence> {
    static void save(diy::BinaryBuffer& bb, const FeatureSequence& m) {
      diy::save(bb, m.its);
      diy::save(bb, m.itl);
      diy::save(bb, m.lids);
      diy::save(bb, m.r);
      diy::save(bb, m.g);
      diy::save(bb, m.b);
    }

    static void load(diy::BinaryBuffer&bb, FeatureSequence& m) {
      diy::load(bb, m.its);
      diy::load(bb, m.itl);
      diy::load(bb, m.lids);
      diy::load(bb, m.r);
      diy::load(bb, m.g);
      diy::load(bb, m.b);
    }
  };
}
#endif

#endif
