#ifndef _MESH_GRAPH_REGULAR_3D_HH
#define _MESH_GRAPH_REGULAR_3D_HH

#include <set>

namespace ftk {

template <typename IndexType>
struct mesh_graph_regular_3d {
  int d[3];
  bool pbc[3];

  // 0, 0: node-node links through edges
  virtual std::set<IndexType> links00(IndexType i) 
  {return std::set<IndexType>();}

  // 0, 1: node-edge links
  virtual std::set<IndexType> links01(IndexType i) 
  {return std::set<IndexType>();}
  
  // 0, 2: node-face links
  virtual std::set<IndexType> links02(IndexType i)
  {return std::set<IndexType>();}
  
  // 0, 3: node-cell links (usually directly available from mesh)
  virtual std::set<IndexType> links03(IndexType i)
  {return std::set<IndexType>();}
  
  // 1, 0: edge-node links
  virtual std::set<IndexType> links10(IndexType i)
  {return std::set<IndexType>();}
  
  // 1, 1: eddg-edge links (an edge links to another edge if they share the same node)
  virtual std::set<IndexType> links11(IndexType i)
  {return std::set<IndexType>();}
  
  // 1, 2: edge-face links
  virtual std::set<IndexType> links12(IndexType i)
  {return std::set<IndexType>();}
  
  // 1, 3: edge-cell links
  virtual std::set<IndexType> links13(IndexType i)
  {return std::set<IndexType>();}
  
  // 2, 0: face-node links
  virtual std::set<IndexType> links20(IndexType i)
  {return std::set<IndexType>();}
  
  // 2, 1: face-edge links
  virtual std::set<IndexType> links21(IndexType i)
  {return std::set<IndexType>();}
  
  // 2, 2: face-face links (a face links to another face if they share the same edge)
  virtual std::set<IndexType> links22(IndexType i)
  {return std::set<IndexType>();}

  // 2, 3: face-cell links
  virtual std::set<IndexType> links23(IndexType i)
  {return std::set<IndexType>();}

  // 3, 0: cell-node links
  virtual std::set<IndexType> links30(IndexType i)
  {return std::set<IndexType>();}

  // 3, 1: cell-edge links
  virtual std::set<IndexType> links31(IndexType i)
  {return std::set<IndexType>();}

  // 3, 2: cell-face links
  virtual std::set<IndexType> links32(IndexType i)
  {return std::set<IndexType>();}

  // 3, 3: cell-cell links (a cell links to another cell if they share the same face)
  virtual std::set<IndexType> links33(IndexType i)
  {return std::set<IndexType>();}
}

}

#endif
