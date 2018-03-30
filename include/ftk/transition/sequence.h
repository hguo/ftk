#ifndef _FTK_COMPONENT_H
#define _FTK_COMPONENT_H

#include <vector>
#include <map>

struct ftkComponent {
  int its, itl;  // start and duration (index of frames)
  std::vector<int> lids; // local ids

  unsigned char r, g, b;

  // std::vector<int> lhs_gids, rhs_gids;
  // int lhs_event, rhs_event;
};

#endif
