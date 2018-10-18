#ifndef _MESHGRAPH_REGULAR3D_H
#define _MESHGRAPH_REGULAR3D_H

#include "common/MeshGraph.h"

class MeshGraphRegular3D : public MeshGraph {
protected:
  int d[3];
  bool pbc[3];

public:
  void nid2nidx(NodeIdType id, int nidx[3]) const;
  NodeIdType nidx2nid(const int nidx[3]) const; // modIdx'ed

  void eid2eidx(EdgeIdType id, int eidx[4]) const;
  EdgeIdType eidx2eid(const int eidx[4]) const;
  
  void fid2fidx(FaceIdType id, int fidx[4]) const;
  FaceIdType fidx2fid(const int fidx[4]) const;

  void cid2cidx(CellIdType id, int cidx[4]) const;
  CellIdType cidx2cid(const int cidx[4]) const;
  
  bool valid_nidx(const int nidx[3]) const;
  bool valid_eidx(const int eidx[4]) const;
  bool valid_fidx(const int fidx[4]) const;
  bool valid_cidx(const int cidx[3]) const;

public:
  MeshGraphRegular3D(int d[3], bool pbc[3]);

  EdgeIdType NEdges() const;
  FaceIdType NFaces() const;
  CellIdType NCells() const;

  CEdge Edge(EdgeIdType i, bool nodes_only=false) const;
  CFace Face(FaceIdType i, bool nodes_only=false) const;
  CCell Cell(CellIdType i, bool nodes_only=false) const;

public:
  std::vector<FaceIdType> GetBoundaryFaceIds(int type) const; // 0: YZ, 1: ZX, 2: XY
};

#endif
