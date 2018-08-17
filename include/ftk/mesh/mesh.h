#ifndef _FTK_MESH
#define _FTK_MESH

#include <set>

namespace ftk {

struct Mesh {
  virtual std::set<size_t> vertexNeighbors(size_t i) const {return std::set<size_t>();}
  virtual std::set<size_t> cellNeighbors(size_t i) const {return std::set<size_t>();}
};

}

#endif
