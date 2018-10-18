#ifndef _MESHGRAPH_REGULAR3DTETS_H
#define _MESHGRAPH_REGULAR3DTETS_H

#include "common/MeshGraphRegular3D.h"

class MeshGraphRegular3DTets : public MeshGraphRegular3D {
public: // nid is the same as regular3D
  void eid2eidx(EdgeIdType id, int eidx[4]) const;
  EdgeIdType eidx2eid(const int eidx[4]) const;
  
  void fid2fidx(FaceIdType id, int fidx[4]) const;
  FaceIdType fidx2fid(const int fidx[4]) const;

  void cid2cidx(CellIdType id, int cidx[4]) const;
  CellIdType cidx2cid(const int cidx[4]) const;
  
  bool valid_eidx(const int eidx[4]) const;
  bool valid_fidx(const int fidx[4]) const;
  bool valid_cidx(const int cidx[4]) const;

public:
  MeshGraphRegular3DTets(int d[3], bool pbc[3]);

  EdgeIdType NEdges() const;
  FaceIdType NFaces() const;
  CellIdType NCells() const;

  CEdge Edge(EdgeIdType i, bool nodes_only=false) const;
  CFace Face(FaceIdType i, bool nodes_only=false) const;
  CCell Cell(CellIdType i, bool nodes_only=false) const;
};

#endif
