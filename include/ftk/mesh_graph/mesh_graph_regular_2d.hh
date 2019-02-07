#ifndef _MESH_GRAPH_REGULAR_2D_HH
#define _MESH_GRAPH_REGULAR_3D_HH

#include <ftk/mesh_graph/mesh_graph.hh>

namespace ftk {

template <typename index_type=size_t, typename chirality_type=signed char>
struct mesh_graph_regular_2d : public mesh_graph<index_type, chirality_type> {
protected: 
  int d[2] = {256, 256}; // dimensions
  bool pbc[2] = {0, 0}; // periodic boundary conditions

public:
  int lb[2] = {0, 0, 0},
      ub[2] = {0, 0, 0};

public:
  mesh_graph_regular_2d(int d_[2]) {
    d[0] = d_[0]; d[1] = d_[1];
    ub[0] = d[0] - 1; ub[1] = d[1] - 1;
  }

public:
  index_type n(int d) const;
  bool valid(int d, index_type i);
  
  void nid2nidx(index_type id, int nidx[2]) const;
  index_type nidx2nid(const int nidx[2]) const; // modIdx'ed

  void eid2eidx(index_type id, int eidx[3]) const;
  index_type eidx2eid(const int eidx[3]) const;
  
  void cid2cidx(index_type id, int cidx[2]) const;
  index_type cidx2cid(const int cidx[2]) const;
  
  bool valid_nidx(const int nidx[2]) const;
  bool valid_eidx(const int eidx[3]) const;
  bool valid_cidx(const int cidx[2]) const;
};

}

#endif
