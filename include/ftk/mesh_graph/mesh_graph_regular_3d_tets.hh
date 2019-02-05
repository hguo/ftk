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
 
  std::vector<std::pair<index_type, chirality_type> > links_face_node(index_type i);
  std::vector<std::pair<index_type, chirality_type> > links_cell_node(index_type i);
  std::vector<std::pair<index_type, chirality_type> > links_node_cell(index_type i);

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

  template <typename T>
  bool locate_point(
      const T x[3], // input coordinates
      index_type nids[4], // output node ids
      T mu[4]); // barycentric coordinates
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
    return 6 * n(0);
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
  fid2fidx(id, idx);
  if (!valid_fidx(idx)) return results;
  
  const int &i = idx[0], &j = idx[1], &k = idx[2], &t = idx[3];
  const int edges_idx[12][3][4] = {
    {{i, j, k, 0}, {i+1, j, k, 2}, {i, j, k, 1}}, 
    {{i, j, k, 1}, {i, j+1, k, 0}, {i, j, k, 2}}, 
    {{i, j, k, 0}, {i+1, j, k, 3}, {i, j, k, 4}}, 
    {{i, j, k, 3}, {i, j, k+1, 0}, {i, j, k, 4}}, 
    {{i, j, k, 2}, {i, j, k, 5}, {i, j, k, 3}}, 
    {{i, j, k, 5}, {i, j, k+1, 2}, {i, j+1, k, 3}}, 
    {{i, j, k, 2}, {i, j, k, 6}, {i, j, k, 4}}, 
    {{i, j, k, 6}, {i+1, j, k+1, 2}, {i, j+1, k, 4}}, 
    {{i, j, k, 5}, {i, j, k+1, 0}, {i, j, k, 6}}, 
    {{i, j+1, k, 0}, {i, j, k, 6}, {i+1, j, k, 5}}, 
    {{i, j, k, 1}, {i+1, j, k, 5}, {i, j, k, 4}}, 
    {{i, j, k, 5}, {i, j, k+1, 1}, {i, j+1, k, 4}}
  }; 
  const chirality_type edges_chi[12][3] = {
    {1, 1, -1}, 
    {1, -1, -1}, 
    {1, 1, -1}, 
    {1, 1, -1}, 
    {1, 1, -1}, 
    {1, 1, -1}, 
    {1, 1, -1}, 
    {1, 1, -1}, 
    {1, 1, -1}, 
    {-1, 1, -1},
    {1, 1, -1}, 
    {1, 1, -1}
  };
  
  for (int p=0; p<3; p++) 
    results.push_back(std::make_pair(eidx2eid(edges_idx[t][p]), edges_chi[t][p]));

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
  const int faces_fidx[6][4][4] = {
    {{i, j, k, 0}, {i, j, k, 2}, {i+1, j, k, 4}, {i, j, k, 10}}, 
    {{i, j, k, 11}, {i, j, k, 5}, {i, j, k+1, 1}, {i, j+1, k, 3}},
    {{i, j, k, 4}, {i, j, k, 6}, {i, j, k, 8}, {i, j, k, 3}}, 
    {{i, j, k, 9}, {i, j+1, k, 2}, {i, j, k, 7}, {i+1, j, k, 5}}, 
    {{i, j, k, 8}, {i, j, k, 11}, {i, j, k+1, 0}, {i, j, k, 7}},
    {{i, j, k, 1}, {i, j, k, 10}, {i, j, k, 9}, {i, j, k, 6}}
  };
  const chirality_type faces_chi[6][4] = {
    {-1, 1, 1, -1}, 
    {-1, 1, 1, 1}, 
    {-1, 1, 1, -1},
    {-1, -1, 1, -1}, 
    {-1, 1, 1, -1}, 
    {-1, 1, 1, -1}
  };
  for (int p=0; p<4; p++)
    results.push_back(std::make_pair(fidx2fid(faces_fidx[t][p]), faces_chi[t][p]));

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
  fid2fidx(id, idx);
  if (!valid_fidx(idx)) return results;

  const int &i = idx[0], &j = idx[1], &k = idx[2], &t = idx[3];
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
  const chirality_type contained_cells_chi[12][2] = {
    {-1, 1}, {-1, 1}, {1, -1}, {-1, 1}, {-1, 1}, {1, -1}, 
    {1, -1}, {-1, 1}, {1, -1}, {-1, 1}, {-1, 1}, {-1, 1}
  };
  
  for (int p=0; p<2; p++) 
    results.push_back(std::make_pair(cidx2cid(contained_cells_cidx[t][p]), contained_cells_chi[t][p]));

  return results;
}

template <typename index_type, typename chirality_type>
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d_tets<index_type, chirality_type>::
links_node_cell(index_type id)
{
  // fprintf(stderr, "calling links_node_cell\n");
  std::vector<std::pair<index_type, chirality_type> > results;
  int nidx[3];
  mesh_graph_regular_3d<index_type, chirality_type>::nid2nidx(id, nidx);
  if (!mesh_graph_regular_3d<index_type, chirality_type>::valid_nidx(nidx)) return results;

  static const int contained_cells[23][4] = {
    {0, 0, 0, 1}, // A, ABCF
    {0, 0, 0, 2}, // A, ADEF
    {0, 0, 0, 5}, // A, ACDF
    {1, 0, 0, 0}, // B, ABCF
    {0, 1, 0, 1}, // D, DEGH
    {0, 1, 0, 2}, // D, ADEF
    {0, 1, 0, 3}, // D, CDFG
    {0, 1, 0, 4}, // D, DEFG
    {0, 1, 0, 5}, // D, ACDF
    {1, 1, 0, 0}, // C, ABCF
    {1, 1, 0, 3}, // C, CDFG
    {1, 1, 0, 5}, // C, ACDF
    {0, 0, 1, 1}, // E, DEGH
    {0, 0, 1, 2}, // E, ADEF
    {0, 0, 1, 4}, // E, DEFG
    {1, 0, 1, 0}, // F, ABCF
    {1, 0, 1, 2}, // F, ADEF
    {1, 0, 1, 3}, // F, CDFG
    {1, 0, 1, 4}, // F, DEFG
    {1, 0, 1, 5}, // F, ACDF
    {1, 1, 1, 1}, // G, DEGH
    {1, 1, 1, 3}, // G, CDFG
    {1, 1, 1, 4}  // G, DEFG
  };

  for (int i = 0; i < 23; i ++) {
    int cidx[4] = {
      -contained_cells[i][0] + nidx[0], 
      -contained_cells[i][1] + nidx[1], 
      -contained_cells[i][2] + nidx[2], 
       contained_cells[i][3]};

    if (valid_cidx(cidx)) {
      index_type cid = cidx2cid(cidx);
      // fprintf(stderr, "adding %d: %d, %d, %d, %d\n", cid, cidx[0], cidx[1], cidx[2], cidx[3]);
      results.push_back(std::make_pair(cid, 0));
    }
  }

  return results;
}

template <typename index_type, typename chirality_type>
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d_tets<index_type, chirality_type>::
links_cell_node(index_type id)
{
  // fprintf(stderr, "calling links_cell_node\n");
  std::vector<std::pair<index_type, chirality_type> > results;
  int cidx[4];
  cid2cidx(id, cidx);
  if (!valid_cidx(cidx)) return results;

  static const int nodes_idx[6][4][3] = {
    {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 0, 1}}, // ABCF
    {{0, 1, 0}, {0, 0, 1}, {1, 1, 1}, {0, 1, 1}}, // DEGH
    {{0, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}}, // ADEF
    {{1, 1, 0}, {0, 1, 0}, {1, 0, 1}, {1, 1, 1}}, // CDFG
    {{0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}}, // DEFG
    {{0, 0, 0}, {1, 1, 0}, {0, 1, 0}, {1, 0, 1}}  // ACDF
  };

  const int type = cidx[3];
  for (int p=0; p<4; p++) {
    int nidx[3];
    for (int q=0; q<3; q++) 
      nidx[q] = cidx[q] + nodes_idx[type][p][q];
    results.push_back(std::make_pair(mesh_graph_regular_3d<index_type, chirality_type>::nidx2nid(nidx), 0)); 
  }

  return results;
}

template <typename index_type, typename chirality_type>
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d_tets<index_type, chirality_type>::
links_face_node(index_type id)
{
  // return mesh_graph<index_type, chirality_type>::links_face_node(id);

  std::vector<std::pair<index_type, chirality_type> > results;
  
  int idx[4];
  fid2fidx(id, idx);
  if (!valid_fidx(idx)) return results;

  const int &i = idx[0], &j = idx[1], &k = idx[2], &t = idx[3];
  static const int nodes_idx[12][3][3] = { // 12 types of faces
    {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}}, // ABC
    {{0, 0, 0}, {1, 1, 0}, {0, 1, 0}}, // ACD
    {{0, 0, 0}, {1, 0, 0}, {1, 0, 1}}, // ABF
    {{0, 0, 0}, {0, 0, 1}, {1, 0, 1}}, // AEF
    {{0, 0, 0}, {0, 1, 0}, {0, 0, 1}}, // ADE
    {{0, 1, 0}, {0, 0, 1}, {0, 1, 1}}, // DEH
    {{0, 0, 0}, {0, 1, 0}, {1, 0, 1}}, // ADF
    {{0, 1, 0}, {1, 0, 1}, {1, 1, 1}}, // DFG
    {{0, 1, 0}, {0, 0, 1}, {1, 0, 1}}, // DEF
    {{1, 1, 0}, {0, 1, 0}, {1, 0, 1}}, // CDF
    {{0, 0, 0}, {1, 1, 0}, {1, 0, 1}}, // ACF
    {{0, 1, 0}, {0, 0, 1}, {1, 1, 1}}  // DEG
  };

  const int type = idx[3];
  for (int p=0; p<3; p++) {
    int nidx[3];
    for (int q=0; q<3; q++) 
      nidx[q] = idx[q] + nodes_idx[type][p][q];
    results.push_back(std::make_pair(mesh_graph_regular_3d<index_type, chirality_type>::nidx2nid(nidx), 0)); 
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
#if 0 // legacy version that takes care of pbc
  static const int nodes_idx[7][2][3] = { // 7 types of edges
    {{0, 0, 0}, {1, 0, 0}}, // AB
    {{0, 0, 0}, {1, 1, 0}}, // AC
    {{0, 0, 0}, {0, 1, 0}}, // AD
    {{0, 0, 0}, {0, 0, 1}}, // AE
    {{0, 0, 0}, {1, 0, 1}}, // AF
    {{0, 1, 0}, {0, 0, 1}}, // DE
    {{0, 1, 0}, {1, 0, 1}}, // DF
  };
  
  for (int dim=0; dim<3; dim++) {
    for (int i=0; i<2; i++) {
      const int p = eidx[dim] + nodes_idx[eidx[3]][i][dim];
      // fprintf(stderr, "type=%d, dim=%d, i=%d, p=%d\n", eidx[3], dim, i, p);
      if (p < 0) return false;
      else if (mesh_graph_regular_3d<index_type, chirality_type>::pbc[dim]) {
        if (p >= mesh_graph_regular_3d<index_type, chirality_type>::d[dim]) 
          return false;
      } else {
        if (p != p % mesh_graph_regular_3d<index_type, chirality_type>::d[dim])
          return false;
      }
    }
  }
  return true;
#else
  const int type = eidx[3]; 
  const int edge_node_offsets[7][2][3] = {
    {{0, 0, 0}, {1, 0, 0}}, 
    {{0, 0, 0}, {1, 1, 0}}, 
    {{0, 0, 0}, {0, 1, 0}}, 
    {{0, 0, 0}, {0, 0, 1}}, 
    {{0, 0, 0}, {1, 0, 1}}, 
    {{0, 1, 0}, {0, 0, 1}}, 
    {{0, 1, 0}, {1, 0, 1}}};

  int nidx0[3], nidx1[3];
  for (int p = 0; p < 3; p ++) {
    nidx0[p] = eidx[p] + edge_node_offsets[type][0][p]; // TODO: pbc
    nidx1[p] = eidx[p] + edge_node_offsets[type][1][p];
  }

  return mesh_graph_regular_3d<index_type, chirality_type>::valid_nidx(nidx0) 
    && mesh_graph_regular_3d<index_type, chirality_type>::valid_nidx(nidx1);
#endif
}

template <typename index_type, typename chirality_type>
bool mesh_graph_regular_3d_tets<index_type, chirality_type>::valid_fidx(const int fidx[4]) const
{
  if (fidx[3]<0 || fidx[3]>=12) return false;
#if 0 // legacy code that handles pbc
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
#else
  const int type = fidx[3];
  static const int face_node_idx[12][3][3] = { // 12 types of faces
    {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}}, // ABC
    {{0, 0, 0}, {1, 1, 0}, {0, 1, 0}}, // ACD
    {{0, 0, 0}, {1, 0, 0}, {1, 0, 1}}, // ABF
    {{0, 0, 0}, {0, 0, 1}, {1, 0, 1}}, // AEF
    {{0, 0, 0}, {0, 1, 0}, {0, 0, 1}}, // ADE
    {{0, 1, 0}, {0, 0, 1}, {0, 1, 1}}, // DEH
    {{0, 0, 0}, {0, 1, 0}, {1, 0, 1}}, // ADF
    {{0, 1, 0}, {1, 0, 1}, {1, 1, 1}}, // DFG
    {{0, 1, 0}, {0, 0, 1}, {1, 0, 1}}, // DEF
    {{1, 1, 0}, {0, 1, 0}, {1, 0, 1}}, // CDF
    {{0, 0, 0}, {1, 1, 0}, {1, 0, 1}}, // ACF
    {{0, 1, 0}, {0, 0, 1}, {1, 1, 1}}  // DEG
  };

  int nidxs[3][3];
  for (int p = 0; p < 3; p ++) {
    for (int q = 0; q < 3; q ++) {
      nidxs[p][q] = fidx[q] + face_node_idx[type][p][q];
    }
  }
  
  return mesh_graph_regular_3d<index_type, chirality_type>::valid_nidx(nidxs[0]) 
    && mesh_graph_regular_3d<index_type, chirality_type>::valid_nidx(nidxs[1])
    && mesh_graph_regular_3d<index_type, chirality_type>::valid_nidx(nidxs[2]);
#endif
}

template <typename index_type, typename chirality_type>
bool mesh_graph_regular_3d_tets<index_type, chirality_type>::valid_cidx(const int idx[3]) const
{
  if (idx[3]<0 || idx[3]>=6) return false; // 6 types of tets

  for (int i=0; i<3; i++)
    if (mesh_graph_regular_3d<index_type, chirality_type>::pbc[i]) {
      if (idx[i] < 0 || idx[i] >= mesh_graph_regular_3d<index_type, chirality_type>::d[i]) return false;
    } else {
      if (idx[i] < mesh_graph_regular_3d<index_type, chirality_type>::lb[i] 
          || idx[i] > mesh_graph_regular_3d<index_type, chirality_type>::ub[i]) return false;
      // if (idx[i] < 0 || idx[i] >= mesh_graph_regular_3d<index_type, chirality_type>::d[i]-1) return false;
    }
  return true;
}

#if 0
template <typename index_type, typename chirality_type>
index_type mesh_graph_regular_3d_tets<index_type, chirality_type>::locate_point(
    const T x[3], index_type nids[4], T mu[4])
{
  int corner[3] = {
    static_cast<int>(x[0]), 
    static_cast<int>(x[1]), 
    static_cast<int>(x[2])
  };

  if (corner[0] < lb[0] || corner[0] > ub[0] 
      || corner[1] < lb[1] || corner[1] > ub[1]
      || corner[2] < lb[2] || corner[2] > ub[2])
    return index_type(-1);

  T y[3] = {
    x[0] - corner[0], 
    x[1] - corner[1], 
    x[2] - corner[2]
  }; 

  // TBA
}
#endif

}

#endif
