#ifndef _MESH_GRAPH_REGULAR_3D_TETS_HH
#define _MESH_GRAPH_REGULAR_3D_TETS_HH

#include <ftk/mesh_graph/mesh_graph_regular_3d.hh>

namespace ftk {

template <typename index_type=size_t, typename chirality_type=signed char>
struct mesh_graph_regular_3d_tets : public mesh_graph_regular_3d<index_type, chirality_type> {
public:
  mesh_graph_regular_3d_tets(int d_[3]) : mesh_graph_regular_3d<index_type, chirality_type>(d_) {}
  mesh_graph_regular_3d_tets(int d_[3], bool pbc_[3]) : mesh_graph_regular_3d<index_type, chirality_type>(d_, pbc_) {}
  mesh_graph_regular_3d_tets(int W, int H, int D) : mesh_graph_regular_3d<index_type, chirality_type>(W, H, D) {}
  mesh_graph_regular_3d_tets(int W, int H, int D, bool pbcx, bool pbcy, bool pbcz) : mesh_graph_regular_3d<index_type, chirality_type>(W, H, D, pbcx, pbcy, pbcz) {}

public:
  index_type n(int d) const;
  
  std::vector<std::pair<index_type, chirality_type> > links_edge_node(index_type i); // returns a list of nodes with chirality=1
  std::vector<std::pair<index_type, chirality_type> > links_face_edge(index_type i); // returns a list of edges with chirality=1
  std::vector<std::pair<index_type, chirality_type> > links_cell_face(index_type i); // returns a list of faces with chirality=1
  
  std::vector<std::pair<index_type, chirality_type> > links_node_edge(index_type i);
  std::vector<std::pair<index_type, chirality_type> > links_edge_face(index_type i);
  std::vector<std::pair<index_type, chirality_type> > links_face_cell(index_type i);
  
  ///////////////////
  bool valid_eidx(const int eidx[4]) const;
  bool valid_fidx(const int fidx[4]) const;
  bool valid_cidx(const int cidx[3]) const;
  
  void eid2eidx(index_type id, int eidx[4]) const;
  index_type eidx2eid(const int eidx[4]) const;
  
  void fid2fidx(index_type id, int fidx[4]) const;
  index_type fidx2fid(const int fidx[4]) const;

  void cid2cidx(index_type id, int cidx[4]) const;
  index_type cidx2cid(const int cidx[4]) const;
};

//////////////////////////////////
template <typename index_type, typename chirality_type>
index_type mesh_graph_regular_3d_tets<index_type, chirality_type>::n(int dim) const
{
  if (dim == 0) { // nodes
    return mesh_graph_regular_3d<index_type, chirality_type>::n(dim);
  } else if (dim == 1) { // edges
    return 7 * n(0); // roughly 7 edges for each node
  } else if (dim == 2) { // faces
    return 12 * n(0); // roughly 12 faces for each node
  } else if (dim == 3) { // cells
    return n(0); // roughly the same number as nodes
  } else return index_type(0);
}

template <typename index_type, typename chirality_type>
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d_tets<index_type, chirality_type>::
links_edge_node(index_type id)
{
  std::vector<std::pair<index_type, chirality_type> > results;
 
  int idx[4];
  eid2eidx(id, idx);
  if (!valid_eidx(idx)) return results;
     
  const int &i = idx[0], &j = idx[1], &k = idx[2], &t = idx[3];
  const int nodes_idx[7][2][3] = {
    {{i, j, k}, {i+1, j, k}},         // AB
    {{i, j, k}, {i+1, j+1, k}},       // AC
    {{i, j, k}, {i, j+1, k}},         // AD
    {{i, j, k}, {i, j, k+1}},         // AE
    {{i, j, k}, {i+1, j, k+1}},       // AF
    {{i, j+1, k}, {i, j, k+1}},       // DE
    {{i, j+1, k}, {i+1, j, k+1}}      // DF
  };
  results.push_back(std::make_pair(mesh_graph_regular_3d<index_type, chirality_type>::nidx2nid(nodes_idx[t][0]), 0)); 
  results.push_back(std::make_pair(mesh_graph_regular_3d<index_type, chirality_type>::nidx2nid(nodes_idx[t][1]), 0));

  return results;
}
  
template <typename index_type, typename chirality_type>
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d_tets<index_type, chirality_type>::
links_face_edge(index_type id)
{
  std::vector<std::pair<index_type, chirality_type> > results;
  
  int idx[4];
  eid2eidx(id, idx);
  if (!valid_eidx(idx)) return results;
  
  const int &i = idx[0], &j = idx[1], &k = idx[2], &t = idx[3];
  const int edges_idx[3][4][4] = {
    {{i, j, k, 1}, {i, j+1, k, 2}, {i, j, k+1, 1}, {i, j, k, 2}},
    {{i, j, k, 2}, {i, j, k+1, 0}, {i+1, j, k, 2}, {i, j, k, 0}},
    {{i, j, k, 0}, {i+1, j, k, 1}, {i, j+1, k, 0}, {i, j, k, 1}}};
  const chirality_type edges_chi[4] = {1, 1, -1, -1};
  for (int p=0; p<4; p++) 
    results.push_back(std::make_pair(eidx2eid(edges_idx[t][p]), edges_chi[p]));

  return results;
}

template <typename index_type, typename chirality_type>
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d_tets<index_type, chirality_type>::
links_cell_face(index_type id)
{
  std::vector<std::pair<index_type, chirality_type> > results;
  
  int idx[4];
  cid2cidx(id, idx);
  if (!valid_cidx(idx)) return results;

  const int &i = idx[0], &j = idx[1], &k = idx[2], &t = idx[3];
  const int faces_fidx[6][4] = {
    {i, j, k, 0}, // type0, yz
    {i, j, k, 1}, // type1, zx
    {i, j, k, 2}, // type2, xy
    {i+1, j, k, 0}, // type0, yz
    {i, j+1, k, 1}, // type1, zx
    {i, j, k+1, 2}};  // type2, xy
  const chirality_type faces_chi[6] = {-1, -1, -1, 1, 1, 1};
  for (int p=0; p<6; p++) 
    results.push_back(std::make_pair(fidx2fid(faces_fidx[p]), faces_chi[p]));

  return results;
}

template <typename index_type, typename chirality_type>
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d_tets<index_type, chirality_type>::
links_node_edge(index_type i)
{
  std::vector<std::pair<index_type, chirality_type> > results;
  // TODO
  return results;
}

template <typename index_type, typename chirality_type>
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d_tets<index_type, chirality_type>::
links_edge_face(index_type id)
{
  std::vector<std::pair<index_type, chirality_type> > results;
  
  int idx[4];
  eid2eidx(id, idx);
  if (!valid_eidx(idx)) return results;
  
  const int &i = idx[0], &j = idx[1], &k = idx[2], &t = idx[3];
  const int contained_faces_fidx[3][4][4] = {
    {{i, j, k, 2}, {i, j, k, 1}, {i, j-1, k, 2}, {i, j, k-1, 1}}, 
    {{i, j, k, 2}, {i, j, k, 0}, {i-1, j, k, 2}, {i, j, k-1, 0}},
    {{i, j, k, 1}, {i, j, k, 0}, {i-1, j, k, 1}, {i, j-1, k, 0}}};
  const chirality_type contained_faces_chi[3][4] = {
    {1, -1, -1, 1}, {-1, 1, 1, -1}, {1, -1, -1, 1}};
  // const int contained_faces_eid[3][4] = {
  //   {0, 3, 2, 1}, {3, 0, 1, 2}, {0, 3, 2, 1}};

  for (int p=0; p<4; p++)
    if (valid_fidx(contained_faces_fidx[t][p]))
      results.push_back(std::make_pair(fidx2fid(contained_faces_fidx[t][p]), contained_faces_chi[t][p]));

  return results;
}

template <typename index_type, typename chirality_type>
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d_tets<index_type, chirality_type>::
links_face_cell(index_type id)
{
  std::vector<std::pair<index_type, chirality_type> > results;
  
  int idx[4];
  cid2cidx(id, idx);
  if (!valid_cidx(idx)) return results;

  const int &i = idx[0], &j = idx[1], &k = idx[2], &t = idx[3];
  const int contained_cells_cidx[3][2][3] = {
    {{i, j, k}, {i-1, j, k}}, 
    {{i, j, k}, {i, j-1, k}},
    {{i, j, k}, {i, j, k-1}}};
  const chirality_type contained_cells_chi[2] = {-1, 1};
  // const int contained_cells_fid[3][2] = {
  //   {0, 3}, {1, 4}, {2, 5}};
  for (int p=0; p<2; p++) {
    results.push_back(std::make_pair(cidx2cid(contained_cells_cidx[t][p]), contained_cells_chi[p]));
    // face.contained_cells.push_back(cidx2cid(contained_cells_cidx[t][p]));
    // face.contained_cells_chirality.push_back(contained_cells_chi[p]);
    // face.contained_cells_fid.push_back(contained_cells_fid[t][p]);
  }

  return results;
}

template <typename index_type, typename chirality_type>
index_type mesh_graph_regular_3d_tets<index_type, chirality_type>::eidx2eid(const int idx[4]) const
{
  return mesh_graph_regular_3d<index_type>::nidx2nid(idx)*7 + idx[3];
}

template <typename index_type, typename chirality_type>
index_type mesh_graph_regular_3d_tets<index_type, chirality_type>::fidx2fid(const int idx[4]) const
{
  return mesh_graph_regular_3d<index_type>::nidx2nid(idx)*12 + idx[3];
}

template <typename index_type, typename chirality_type>
index_type mesh_graph_regular_3d_tets<index_type, chirality_type>::cidx2cid(const int idx[3]) const
{
  return mesh_graph_regular_3d<index_type>::nidx2nid(idx)*6 + idx[3];
}

template <typename index_type, typename chirality_type>
void mesh_graph_regular_3d_tets<index_type, chirality_type>::eid2eidx(index_type id, int idx[4]) const
{
  unsigned int nid = id / 7;
  mesh_graph_regular_3d<index_type>::nid2nidx(nid, idx);
  idx[3] = id % 7;
}

template <typename index_type, typename chirality_type>
void mesh_graph_regular_3d_tets<index_type, chirality_type>::fid2fidx(index_type id, int idx[4]) const
{
  unsigned int nid = id / 12;
  mesh_graph_regular_3d<index_type>::nid2nidx(nid, idx);
  idx[3] = id % 12;
}

template <typename index_type, typename chirality_type>
void mesh_graph_regular_3d_tets<index_type, chirality_type>::cid2cidx(index_type id, int idx[3]) const
{
  unsigned int nid = id / 6;
  mesh_graph_regular_3d<index_type>::nid2nidx(nid, idx);
  idx[3] = id % 6;
}

template <typename index_type, typename chirality_type>
bool mesh_graph_regular_3d_tets<index_type, chirality_type>::valid_eidx(const int eidx[4]) const
{
  if (eidx[3]<0 || eidx[3]>=7) return false;
  else {
    for (int i=0; i<3; i++)
      if (mesh_graph_regular_3d<index_type, chirality_type>::pbc[i]) {
        if (eidx[i] < 0 || eidx[i] >= mesh_graph_regular_3d<index_type, chirality_type>::d[i]) return false;
      } else {
        if (eidx[i] < 0 || eidx[i] >= mesh_graph_regular_3d<index_type, chirality_type>::d[i]-1) return false;
      }
    return true;
  }
}

template <typename index_type, typename chirality_type>
bool mesh_graph_regular_3d_tets<index_type, chirality_type>::valid_fidx(const int fidx[4]) const
{
  if (fidx[3]<0 || fidx[3]>=12) return false;
  else {
    int o[3] = {0};
    for (int i=0; i<3; i++)
      if (mesh_graph_regular_3d<index_type, chirality_type>::pbc[i]) {
        if (fidx[i] < 0 || fidx[i] >= mesh_graph_regular_3d<index_type, chirality_type>::d[i]) return false;
      } else {
        if (fidx[i] < 0 || fidx[i] > mesh_graph_regular_3d<index_type, chirality_type>::d[i]-1) return false;
        else if (fidx[i] == mesh_graph_regular_3d<index_type, chirality_type>::d[i]-1) o[i] = 1;
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

template <typename index_type, typename chirality_type>
bool mesh_graph_regular_3d_tets<index_type, chirality_type>::valid_cidx(const int idx[3]) const
{
  if (idx[3]<0 || idx[3]>=6) return false; // 6 types of tets

  for (int i=0; i<3; i++)
    if (mesh_graph_regular_3d<index_type, chirality_type>::pbc[i]) {
      if (idx[i] < 0 || idx[i] >= mesh_graph_regular_3d<index_type, chirality_type>::d[i]) return false;
    } else {
      if (idx[i] < 0 || idx[i] >= mesh_graph_regular_3d<index_type, chirality_type>::d[i]-1) return false;
    }
  return true;
}

}

#endif
