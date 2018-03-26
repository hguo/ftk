#ifndef _FTK_SEQUENCE_H
#define _FTK_SEQUENCE_H

#include <vector>
#include <map>
#include "common/ftkTransition.h"

struct ftkSequence {
  // int ts, tl; // start and duration
  int its, itl;  // start and duration (index of frames)
  std::vector<int> lids; // local ids

  unsigned char r, g, b;

  // std::vector<int> lhs_gids, rhs_gids;
  // int lhs_event, rhs_event;
};

#endif
