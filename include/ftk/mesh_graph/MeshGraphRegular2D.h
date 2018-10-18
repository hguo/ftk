#ifndef _MESHGRAPH_REGULAR2D_H
#define _MESHGRAPH_REGULAR2D_H

#include "common/MeshGraph.h"

class MeshGraphRegular2D : public MeshGraph {
private:
  int d[2];
  bool pbc[2];

private:
  void nid2nidx(NodeIdType id, int nidx[2]) const;
  NodeIdType nidx2nid(const int nidx[2]) const; // modIdx'ed

  void eid2eidx(EdgeIdType id, int eidx[3]) const;
  EdgeIdType eidx2eid(const int eidx[3]) const;
  
  void fid2fidx(FaceIdType id, int fidx[2]) const;
  FaceIdType fidx2fid(const int fidx[2]) const;

  void cid2cidx(CellIdType id, int cidx[2]) const;
  CellIdType cidx2cid(const int cidx[2]) const;
  
  bool valid_nidx(const int nidx[2]) const;
  bool valid_eidx(const int eidx[3]) const;
  bool valid_fidx(const int fidx[3]) const;
  bool valid_cidx(const int cidx[2]) const;

public:
  MeshGraphRegular2D(int d[3], bool pbc[3]);
  
  EdgeIdType NEdges() const;
  FaceIdType NFaces() const;
  CellIdType NCells() const;

  CEdge Edge(EdgeIdType i, bool nodes_only=false) const;
  CFace Face(FaceIdType i, bool nodes_only=false) const;
  CCell Cell(CellIdType i, bool nodes_only=false) const;
};

#endif
