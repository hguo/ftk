#include "MeshGraphRegular2D.h"
#include <cstring>
#include <climits>

MeshGraphRegular2D::MeshGraphRegular2D(int d_[2], bool pbc_[2])
{
  memcpy(d, d_, sizeof(int)*2);
  memcpy(pbc, pbc_, sizeof(bool)*2);
  // pbc[0] = pbc[1] = 0; // disable pbc in the test
}

CCell MeshGraphRegular2D::Cell(CellIdType id, bool nodes_only) const
{
  CCell cell;
  int idx[2];
  
  cid2cidx(id, idx);
  if (!valid_cidx(idx)) return cell;
  const int i = idx[0], j = idx[1];

  // nodes
  const int nodes_idx[4][2] = {{i, j}, {i+1, j}, {i+1, j+1}, {i, j+1}};
  for (int p=0; p<4; p++) // don't worry about modIdx here. automatically done in idx2id()
    cell.nodes.push_back(nidx2nid(nodes_idx[i]));
  if (nodes_only) return cell;

  // faces
  const int face_fidx[2] = {i, j}; // only one face in the cell
  const ChiralityType face_chi = 1;
  cell.faces.push_back(fidx2fid(face_fidx));
  cell.faces_chirality.push_back(face_chi);

  // neighbor cells (no neighbor cells for the face)
  cell.neighbor_cells.push_back(UINT_MAX);

  return cell;
}

CFace MeshGraphRegular2D::Face(FaceIdType id, bool nodes_only) const
{
  CFace face;
  int fidx[3];

  fid2fidx(id, fidx);
  if (!valid_fidx(fidx)) return face;
  const int i = fidx[0], j = fidx[1];

  // nodes
  const int nodes_idx[4][2] = {{i, j}, {i+1, j}, {i+1, j+1}, {i, j+1}};
  for (int p=0; p<4; p++) // don't worry about modIdx here. automatically done in idx2id()
    face.nodes.push_back(nidx2nid(nodes_idx[p]));
  if (nodes_only) return face;

  // edges
  const int edges_idx[4][3] = {{i, j, 0}, {i+1, j, 1}, {i, j+1, 0}, {i, j, 1}};
  const ChiralityType edges_chi[4] = {1, 1, -1, -1};
  for (int p=0; p<4; p++) {
    face.edges.push_back(eidx2eid(edges_idx[p]));
    face.edges_chirality.push_back(edges_chi[p]);
  }

  // contained cells
  const int contained_cells_cidx[2] = {i, j};
  const ChiralityType contained_cells_chi = 1;
  const int contained_cells_fid = 0;
  face.contained_cells.push_back(cidx2cid(contained_cells_cidx));
  face.contained_cells_chirality.push_back(contained_cells_chi);
  face.contained_cells_fid.push_back(contained_cells_fid);

  return face;
}

CEdge MeshGraphRegular2D::Edge(EdgeIdType id, bool nodes_only) const
{
  CEdge edge;
  int eidx[3];

  eid2eidx(id, eidx);
  if (!valid_eidx(eidx)) return edge;
  const int i = eidx[0], j = eidx[1], t = eidx[2];

  // nodes
  const int nodes_idx[2][2][2] = {
    {{i, j}, {i+1, j}}, 
    {{i, j}, {i, j+1}}};

  edge.node0 = nidx2nid(nodes_idx[t][0]); // id; // nidx2nid(node_idx[t][0]);
  edge.node1 = nidx2nid(nodes_idx[t][1]);
  if (nodes_only) return edge;

  // contained faces
  const int contained_faces_fidx[2][2][2] = {
    {{i, j}, {i, j-1}}, 
    {{i, j}, {i-1, j}}};
  const ChiralityType contained_faces_chi[2][2] = {
    {1, -1}, {-1, 1}};
  const int contained_faces_eid[2][2] = {
    {0, 2}, {3, 1}};
  for (int p=0; p<2; p++) {
    edge.contained_faces.push_back(fidx2fid(contained_faces_fidx[t][p]));
    edge.contained_faces_chirality.push_back(contained_faces_chi[t][p]);
    edge.contained_faces_eid.push_back(contained_faces_eid[t][p]);
  }

  return edge;
}

EdgeIdType MeshGraphRegular2D::NEdges() const
{
  return NCells()*2;
}

EdgeIdType MeshGraphRegular2D::NFaces() const
{
  return NCells();
}

EdgeIdType MeshGraphRegular2D::NCells() const
{
  return d[0]*d[1];
}

void MeshGraphRegular2D::nid2nidx(unsigned int id, int idx[2]) const
{
  int j = id / d[0]; 
  int i = id - j*d[0];

  idx[0] = i; idx[1] = j;
}

void MeshGraphRegular2D::eid2eidx(unsigned int id, int idx[3]) const
{
  const unsigned int nid = id / 2;
  nid2nidx(nid, idx);
  idx[2] = id % 2;
}

void MeshGraphRegular2D::fid2fidx(unsigned int id, int idx[2]) const
{
  unsigned int nid = id;
  nid2nidx(id, idx);
}

void MeshGraphRegular2D::cid2cidx(unsigned int id, int idx[2]) const
{
  nid2nidx(id, idx);
}

unsigned int MeshGraphRegular2D::nidx2nid(const int idx_[2]) const
{
  int idx[2] = {idx_[0], idx_[1]};
  for (int i=0; i<2; i++) {
    idx[i] = idx[i] % d[i];
    if (idx[i] < 0)
      idx[i] += d[i];
  }
  return idx[0] + d[0] * idx[1];
}

unsigned int MeshGraphRegular2D::eidx2eid(const int idx[3]) const
{
  return nidx2nid(idx)*2 + idx[2];
}

unsigned int MeshGraphRegular2D::fidx2fid(const int idx[3]) const
{
  return nidx2nid(idx);
}

unsigned int MeshGraphRegular2D::cidx2cid(const int idx[2]) const
{
  return nidx2nid(idx);
}

bool MeshGraphRegular2D::valid_nidx(const int idx[2]) const
{
  for (int i=0; i<2; i++) 
    if (idx[i]<0 || idx[i]>=d[i]) 
      return false;
  return true;
}

bool MeshGraphRegular2D::valid_eidx(const int eidx[3]) const
{
  if (eidx[2]<0 || eidx[2]>=2) return false;
  else return valid_cidx(eidx);
}

bool MeshGraphRegular2D::valid_fidx(const int fidx[3]) const
{
  return valid_cidx(fidx);
}

bool MeshGraphRegular2D::valid_cidx(const int idx[2]) const
{
  for (int i=0; i<2; i++)
    if (pbc[i]) {
      if (idx[i] >= d[i]) return false;
    } else {
      if (idx[i] >= d[i]-1) return false;
    }
  return true;
}
