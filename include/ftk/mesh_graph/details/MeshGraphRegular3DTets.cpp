#include "MeshGraphRegular3DTets.h"
#include <cstring>
#include <climits>

MeshGraphRegular3DTets::MeshGraphRegular3DTets(int d[3], bool pbc[3])
  : MeshGraphRegular3D(d, pbc)
{
}

CCell MeshGraphRegular3DTets::Cell(CellIdType id, bool nodes_only) const
{
  CCell cell; 
  int idx[4];

  cid2cidx(id, idx);
  if (!valid_cidx(idx)) return cell;
  const int i = idx[0], j = idx[1], k = idx[2], t = idx[3];

  const int nodes_idx[6][4][3] = {
    {{i, j, k}, {i+1, j, k}, {i+1, j+1, k}, {i+1, j, k+1}}, 
    {{i, j+1, k}, {i, j, k+1}, {i+1, j+1, k+1}, {i, j+1, k+1}}, 
    {{i, j, k}, {i, j+1, k}, {i, j, k+1}, {i+1, j, k+1}}, 
    {{i+1, j+1, k}, {i, j+1, k}, {i+1, j, k+1}, {i+1, j+1, k+1}}, 
    {{i, j+1, k}, {i, j, k+1}, {i+1, j, k+1}, {i+1, j+1, k+1}}, 
    {{i, j, k}, {i+1, j+1, k}, {i, j+1, k}, {i+1, j, k+1}}
  };
  
  for (int p=0; p<4; p++) 
    cell.nodes.push_back(nidx2nid(nodes_idx[t][p]));
  if (nodes_only) return cell;

  // faces
  const int faces_fidx[6][4][4] = {
    {{i, j, k, 0}, {i, j, k, 2}, {i+1, j, k, 4}, {i, j, k, 10}}, 
    {{i, j, k, 11}, {i, j, k, 5}, {i, j, k+1, 1}, {i, j+1, k, 3}},
    {{i, j, k, 4}, {i, j, k, 6}, {i, j, k, 8}, {i, j, k, 3}}, 
    {{i, j, k, 9}, {i, j+1, k, 2}, {i, j, k, 7}, {i+1, j, k, 5}}, 
    {{i, j, k, 8}, {i, j, k, 11}, {i, j, k+1, 0}, {i, j, k, 7}},
    {{i, j, k, 1}, {i, j, k, 10}, {i, j, k, 9}, {i, j, k, 6}}
  };
  const ChiralityType faces_chi[6][4] = {
    {-1, 1, 1, -1}, 
    {-1, 1, 1, 1}, 
    {-1, 1, 1, -1},
    {-1, -1, 1, -1}, 
    {-1, 1, 1, -1}, 
    {-1, 1, 1, -1}
  };
  for (int p=0; p<4; p++) {
    cell.faces.push_back(fidx2fid(faces_fidx[t][p]));
    cell.faces_chirality.push_back(faces_chi[t][p]);
  }
  
  // neighbor cells
  const int neighbors_cidx[6][4][4] = { // need to be consistent with faces
    {{i, j, k-1, 4}, {i, j-1, k, 3}, {i+1, j, k, 2}, {i, j, k, 5}}, 
    {{i, j, k, 4}, {i-1, j, k, 3}, {i, j, k+1, 5}, {i, j+1, k, 2}}, 
    {{i-1, j, k, 0}, {i, j, k, 5}, {i, j, k, 4}, {i, j-1, k, 1}}, 
    {{i, j, k, 5}, {i, j+1, k, 0}, {i, j, k, 4}, {i+1, j, k, 1}}, 
    {{i, j, k, 2}, {i, j, k, 1}, {i, j, k+1, 0}, {i, j, k, 3}},
    {{i, j, k-1, 1}, {i, j, k, 0}, {i, j, k, 3}, {i, j, k, 2}}
  }; 
  for (int p=0; p<4; p++)
    cell.neighbor_cells.push_back(cidx2cid(neighbors_cidx[t][p]));

  return cell;
}

CFace MeshGraphRegular3DTets::Face(FaceIdType id, bool nodes_only) const
{
  CFace face;
  int fidx[4];

  fid2fidx(id, fidx);
  bool valid = valid_fidx(fidx);
  const int i = fidx[0], j = fidx[1], k = fidx[2], t = fidx[3];
  if (!valid) {
    // fprintf(stderr, "invalid {%d, %d, %d, %d}\n", i, j, k, t);
    return face;
  }

  // nodes
  const int nodes_idx[12][3][3] = { // 12 types of faces
    {{i, j, k}, {i+1, j, k}, {i+1, j+1, k}},      // 0: ABC
    {{i, j, k}, {i+1, j+1, k}, {i, j+1, k}},      // 1: ACD
    {{i, j, k}, {i+1, j, k}, {i+1, j, k+1}},      // 2: ABF
    {{i, j, k}, {i, j, k+1}, {i+1, j, k+1}},      // 3: AEF
    {{i, j, k}, {i, j+1, k}, {i, j, k+1}},        // 4: ADE
    {{i, j+1, k}, {i, j, k+1}, {i, j+1, k+1}},    // 5: DEH
    {{i, j, k}, {i, j+1, k}, {i+1, j, k+1}},      // 6: ADF
    {{i, j+1, k}, {i+1, j, k+1}, {i+1, j+1, k+1}},// 8: DFG
    {{i, j+1, k}, {i, j, k+1}, {i+1, j, k+1}},    // 9: DEF
    {{i+1, j+1, k}, {i, j+1, k}, {i+1, j, k+1}},  //10: CDF
    {{i, j, k}, {i+1, j+1, k}, {i+1, j, k+1}},    //11: ACF
    {{i, j+1, k}, {i, j, k+1}, {i+1, j+1, k+1}}   //12: DEG
  };
  for (int p=0; p<3; p++)
    face.nodes.push_back(nidx2nid(nodes_idx[t][p]));
  if (nodes_only) return face;

  // edges
  const int edges_idx[12][3][4] = {
    {{i, j, k, 0}, {i+1, j, k, 2}, {i, j, k, 1}}, 
    {{i, j, k, 1}, {i, j+1, k, 0}, {i, j, k, 2}}, 
    {{i, j, k, 0}, {i+1, j, k, 3}, {i, j, k, 4}}, 
    {{i, j, k, 3}, {i, j, k+1, 0}, {i, j, k, 4}}, 
    {{i, j, k, 2}, {i, j, k, 5}, {i, j, k, 3}}, 
    {{i, j, k, 5}, {i, j, k, 2}, {i, j+1, k, 3}}, 
    {{i, j, k, 2}, {i, j, k, 6}, {i, j, k, 4}}, 
    {{i, j, k, 6}, {i+1, j, k+1, 2}, {i, j+1, k, 4}}, 
    {{i, j, k, 5}, {i, j, k+1, 0}, {i, j, k, 6}}, 
    {{i, j+1, k, 0}, {i, j, k, 6}, {i+1, j, k, 5}}, 
    {{i, j, k, 1}, {i+1, j, k, 5}, {i, j, k, 4}}, 
    {{i, j, k, 5}, {i, j, k+1, 1}, {i, j+1, k, 4}}
  }; 
  const ChiralityType edges_chi[12][3] = {
    {1, 1, -1}, 
    {1, -1, -1}, 
    {1, 1, -1}, 
    {1, 1, -1}, 
    {1, 1, -1}, 
    {1, 1, -1}, 
    {1, 1, -1}, 
    {1, 1, -1}, 
    {-1, 1, -1},
    {1, -1, -1}, 
    {1, 1, -1}
  };
    
  for (int p=0; p<3; p++) {
    face.edges.push_back(eidx2eid(edges_idx[t][p]));
    face.edges_chirality.push_back(edges_chi[t][p]);
  }

  // contained cells
  const int contained_cells_cidx[12][2][4] = {
    {{i, j, k, 0}, {i, j, k-1, 4}},  // ABC
    {{i, j, k, 5}, {i, j, k-1, 1}},  // ACD
    {{i, j, k, 0}, {i, j-1, k, 3}},  // ABF
    {{i, j, k, 2}, {i, j-1, k, 1}},  // AEF
    {{i, j, k, 2}, {i-1, j, k, 0}},  // ADE
    {{i, j, k, 1}, {i-1, j, k, 3}},  // DEH
    {{i, j, k, 2}, {i, j, k, 5}},    // ADF
    {{i, j, k, 4}, {i, j, k, 3}},    // DFG
    {{i, j, k, 2}, {i, j, k, 4}},    // DEF
    {{i, j, k, 3}, {i, j, k, 5}},    // CDF
    {{i, j, k, 0}, {i, j, k, 5}},    // ACF
    {{i, j, k, 1}, {i, j, k, 4}}     // DEG
  }; 
  const ChiralityType contained_cells_chi[12][2] = {
    {-1, 1}, 
    {-1, 1}, 
    {1, -1}, 
    {-1, 1}, 
    {-1, 1}, 
    {1, -1},
    {1, -1}, 
    {-1, 1}, 
    {1, -1}, 
    {-1, 1}, 
    {-1, 1}, 
    {-1, 1}
  };
  const int contained_cells_fid[12][2] = {
    {0, 2}, 
    {0, 2}, 
    {1, 1},
    {3, 3}, 
    {0, 2}, 
    {1, 3}, 
    {1, 3}, 
    {3, 2}, 
    {2, 0}, 
    {0, 2},
    {3, 1}, 
    {0, 1}
  };
  for (int p=0; p<2; p++) {
    if (!valid_cidx(contained_cells_cidx[t][p])) continue;
    face.contained_cells.push_back(cidx2cid(contained_cells_cidx[t][p]));
    face.contained_cells_chirality.push_back(contained_cells_chi[t][p]);
    face.contained_cells_fid.push_back(contained_cells_fid[t][p]);
  }

#if 0
  fprintf(stderr, "fid=%u, fidx={%d, %d, %d, %d}, cid0=%u={%d, %d, %d, %d}, fid0=%d, cid1=%u={%d, %d, %d, %d}, fid1=%d\n", 
      id, i, j, k, t, 
      face.contained_cells[0], contained_cells_cidx[t][0][0], contained_cells_cidx[t][0][1], contained_cells_cidx[t][0][2], contained_cells_cidx[t][0][3], contained_cells_fid[t][0], 
      face.contained_cells[1], contained_cells_cidx[t][1][0], contained_cells_cidx[t][1][1], contained_cells_cidx[t][1][2], contained_cells_cidx[t][1][3], contained_cells_fid[t][1]);
#endif

  return face;
}

CEdge MeshGraphRegular3DTets::Edge(EdgeIdType id, bool nodes_only) const
{
  CEdge edge;
  int eidx[4];

  eid2eidx(id, eidx);
  if (!valid_eidx(eidx)) return edge;
  const int i = eidx[0], j = eidx[1], k = eidx[2], t = eidx[3];

  // nodes
  const int nodes_idx[7][2][3] = {
    {{i, j, k}, {i+1, j, k}},         // AB
    {{i, j, k}, {i+1, j+1, k}},       // AC
    {{i, j, k}, {i, j+1, k}},         // AD
    {{i, j, k}, {i, j, k+1}},         // AE
    {{i, j, k}, {i+1, j, k+1}},       // AF
    {{i, j+1, k}, {i, j, k+1}},       // DE
    {{i, j+1, k}, {i+1, j, k+1}}      // DF
  };

  edge.node0 = nidx2nid(nodes_idx[t][0]);
  edge.node1 = nidx2nid(nodes_idx[t][1]);
  if (nodes_only) return edge;

  // contained faces (each edge connects to 4 or 6 faces)
  const int contained_faces_fidx[7][6][4] = {
    {{i, j, k, 0}, {i, j, k, 2}, {i, j-1, k, 9}, {i, j-1, k, 1}, {i, j, k-1, 8}, {i, j, k-1, 3}}, 
    {{i, j, k, 0}, {i, j, k, 1}, {i, j, k, 10}, {i, j, k-1, 11}, {0, 0, 0, -1}, {0, 0, 0, -1}}, 
    {{i, j, k, 1}, {i, j, k, 6}, {i, j, k, 4}, {i-1, j, k, 0}, {i, j, k-1, 5}, {i-1, j, k-1, 7}}, 
    {{i, j, k, 3}, {i, j, k, 4}, {i-1, j, k, 2}, {i, j-1, k, 5}, {0, 0, 0, -1}, {0, 0, 0, -1}}, 
    {{i, j, k, 2}, {i, j, k, 10}, {i, j, k, 6}, {i, j, k, 3}, {i, j-1, k, 11}, {i, j-1, k, 7}}, 
    {{i, j, k, 8}, {i, j, k, 11}, {i, j, k, 5}, {i, j, k, 4}, {i-1, j, k, 10}, {i-1, j, k, 9}}, 
    {{i, j, k, 6}, {i, j, k, 8}, {i, j, k, 7}, {i, j, k, 9}, {0, 0, 0, -1}, {0, 0, 0, -1}}
  };
  const ChiralityType contained_faces_chi[7][6] = {
    {1, 1, -1, -1, 1, 1}, 
    {-1, 1, 1, 1, 0, 0}, 
    {1, 1, 1, 1, 1, 1},
    {1, -1, 1, -1, 0, 0}, 
    {-1, -1, -1, -1, -1, -1}, 
    {1, 1, 1, 1, 1, -1}, 
    {1, -1, 1, 1, 0, 0}
  };
  const int contained_faces_eid[7][6] = {
    {0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 1, -1, -1}, 
    {2, 0, 0, 1, 1, 1}, 
    {0, 2, 1, 2, -1, -1}, 
    {2, 2, 2, 2, 2, 2}, 
    {0, 0, 0, 1, 1, 2}, 
    {1, 2, 0, 1, -1, -1}
  }; 

  for (int p=0; p<6; p++) {
    if (contained_faces_chi[t][p] != 0) {
      edge.contained_faces.push_back(fidx2fid(contained_faces_fidx[t][p]));
      edge.contained_faces_chirality.push_back(contained_faces_chi[t][p]);
      edge.contained_faces_eid.push_back(contained_faces_eid[t][p]);
    }
  }

  return edge;
}

EdgeIdType MeshGraphRegular3DTets::NEdges() const
{
  return NCells()*7;
}

EdgeIdType MeshGraphRegular3DTets::NFaces() const
{
  return NCells()*12;
}

EdgeIdType MeshGraphRegular3DTets::NCells() const
{
  return d[0]*d[1]*d[2];
}

void MeshGraphRegular3DTets::eid2eidx(unsigned int id, int idx[4]) const
{
  unsigned int nid = id / 7;
  nid2nidx(nid, idx);
  idx[3] = id % 7;
}

void MeshGraphRegular3DTets::fid2fidx(unsigned int id, int idx[4]) const
{
  unsigned int nid = id / 12;
  nid2nidx(nid, idx);
  idx[3] = id % 12;
}

void MeshGraphRegular3DTets::cid2cidx(unsigned int id, int idx[4]) const
{
  unsigned int nid = id / 6;
  nid2nidx(nid, idx);
  idx[3] = id % 6;
}

unsigned int MeshGraphRegular3DTets::eidx2eid(const int idx[4]) const
{
  return nidx2nid(idx)*7 + idx[3];
}

unsigned int MeshGraphRegular3DTets::fidx2fid(const int idx[4]) const
{
  return nidx2nid(idx)*12 + idx[3];
}

unsigned int MeshGraphRegular3DTets::cidx2cid(const int idx[4]) const
{
  return nidx2nid(idx)*6 + idx[3];
}

bool MeshGraphRegular3DTets::valid_eidx(const int eidx[4]) const
{
  if (eidx[3]<0 || eidx[3]>=7) return false;
  else {
    for (int i=0; i<3; i++)
      if (pbc[i]) {
        if (eidx[i] < 0 || eidx[i] >= d[i]) return false;
      } else {
        if (eidx[i] < 0 || eidx[i] >= d[i]-1) return false;
      }
    return true;
  }
}

bool MeshGraphRegular3DTets::valid_fidx(const int fidx[4]) const
{
  if (fidx[3]<0 || fidx[3]>=12) return false;
  else {
    int o[3] = {0};
    for (int i=0; i<3; i++)
      if (pbc[i]) {
        if (fidx[i] < 0 || fidx[i] >= d[i]) return false;
      } else {
        if (fidx[i] < 0 || fidx[i] > d[i]-1) return false;
        else if (fidx[i] == d[i]-1) o[i] = 1;
      }
    
    const int sum = o[0] + o[1] + o[2];
    if (sum == 0) return true;
    else if (o[0] + o[1] + o[2] > 1) return false;
    else if (o[0] && (fidx[3] == 4 || fidx[3] == 5)) return true;
    else if (o[1] && (fidx[3] == 2 || fidx[3] == 3)) return true; 
    else if (o[2] && (fidx[3] == 0 || fidx[3] == 1)) return true;
    else return false;
  }
}

bool MeshGraphRegular3DTets::valid_cidx(const int idx[4]) const
{
  if (idx[3]<0 || idx[3]>=6) return false; // 6 types of tets

  for (int i=0; i<3; i++)
    if (pbc[i]) {
      if (idx[i] < 0 || idx[i] >= d[i]) return false;
    } else {
      if (idx[i] < 0 || idx[i] >= d[i]-1) return false;
    }
  return true;
}

