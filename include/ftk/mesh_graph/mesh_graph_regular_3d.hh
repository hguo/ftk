#ifndef _MESH_GRAPH_REGULAR_3D_HH
#define _MESH_GRAPH_REGULAR_3D_HH

// #include <set>
#include <vector>
#include <ftk/mesh_graph/mesh_graph.hh>

namespace ftk {

template <typename index_type=size_t, typename chirality_type=signed char>
struct mesh_graph_regular_3d : public mesh_graph<index_type, chirality_type> {
protected:
  int d[3] = {256, 256, 256}; // dimensions
  bool pbc[3] = {0, 0, 0}; // periodic boundary conditions

public:
  mesh_graph_regular_3d(int d_[3]) {
    d[0] = d_[0]; d[1] = d_[1]; d[2] = d_[2];
  }
  
  mesh_graph_regular_3d(int d_[3], bool pbc_[3]) {
    d[0] = d_[0]; d[1] = d_[1]; d[2] = d_[2];
    pbc[0] = pbc_[0]; pbc[1] = pbc_[1]; pbc[2] = pbc_[2];
  }

  mesh_graph_regular_3d(int W, int H, int D) {
    d[0] = W; d[1] = H; d[2] = D;
  }

  mesh_graph_regular_3d(int W, int H, int D, bool pbcx, bool pbcy, bool pbcz) :
    mesh_graph_regular_3d(W, H, D) {
    pbc[0] = pbcx; pbc[1] = pbcy; pbc[2] = pbcz;
  }

public:
  index_type n(int d) const;
  bool valid(int d, index_type i);
  
  std::vector<std::pair<index_type, chirality_type> > links_edge_node(index_type i); // returns a list of nodes with chirality=1
  std::vector<std::pair<index_type, chirality_type> > links_face_edge(index_type i); // returns a list of edges with chirality=1
  std::vector<std::pair<index_type, chirality_type> > links_cell_face(index_type i); // returns a list of faces with chirality=1
  
  std::vector<std::pair<index_type, chirality_type> > links_node_edge(index_type i);
  std::vector<std::pair<index_type, chirality_type> > links_edge_face(index_type i);
  std::vector<std::pair<index_type, chirality_type> > links_face_cell(index_type i);
 
  ///////////////////
  bool valid_nidx(const int nidx[3]) const;
  virtual bool valid_eidx(const int eidx[4]) const;
  virtual bool valid_fidx(const int fidx[4]) const;
  virtual bool valid_cidx(const int cidx[3]) const;
  
  virtual void nid2nidx(index_type id, int nidx[3]) const;
  virtual index_type nidx2nid(const int nidx[3]) const; // modIdx'ed

  virtual void eid2eidx(index_type id, int eidx[4]) const;
  virtual index_type eidx2eid(const int eidx[4]) const;
  
  virtual void fid2fidx(index_type id, int fidx[4]) const;
  virtual index_type fidx2fid(const int fidx[4]) const;

  virtual void cid2cidx(index_type id, int cidx[4]) const;
  virtual index_type cidx2cid(const int cidx[4]) const;
};

////////////////////////////
template <typename index_type, typename chirality_type>
index_type mesh_graph_regular_3d<index_type, chirality_type>::n(int dim) const
{
  if (dim == 0) { // nodes
    return index_type(d[0]) * index_type(d[1]) * index_type(d[2]);
  } else if (dim == 1) { // edges
    return 3 * n(0); // roughly 3 edges for each node
  } else if (dim == 2) { // faces
    return 3 * n(0); // roughly 3 faces for each node
  } else if (dim == 3) { // cells
    return n(0); // roughly the same number as nodes
  } else return index_type(0);
}

template <typename index_type, typename chirality_type>
bool mesh_graph_regular_3d<index_type, chirality_type>::valid(int dim, index_type i)
{
  int idx[4];
  if (dim == 0) {
    nid2nidx(i, idx);
    return valid_nidx(idx);
  } else if (dim == 1) {
    eid2eidx(i, idx);
    return valid_eidx(idx);
  } else if (dim == 2) {
    fid2fidx(i, idx);
    return valid_fidx(idx);
  } else if (dim == 3) {
    cid2cidx(i, idx);
    return valid_cidx(idx);
  } else return false;
}

template <typename index_type, typename chirality_type>
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d<index_type, chirality_type>::
links_edge_node(index_type id)
{
  std::vector<std::pair<index_type, chirality_type> > results;

  int idx[4];
  eid2eidx(id, idx);
  // fprintf(stderr, "-eid=%lu, eidx={%d, %d, %d, %d}, valid=%d\n",
  //     id, idx[0], idx[1], idx[2], idx[3], valid_eidx(idx));

  if (!valid_eidx(idx)) return results;
     
  const int &i = idx[0], &j = idx[1], &k = idx[2], &type = idx[3];
  const int nodes_idx[3][2][3] = {
    {{i, j, k}, {i+1, j, k}}, 
    {{i, j, k}, {i, j+1, k}}, 
    {{i, j, k}, {i, j, k+1}}
  };

  results.push_back(std::make_pair(nidx2nid(nodes_idx[type][0]), 0)); // node chirality is 0
  results.push_back(std::make_pair(nidx2nid(nodes_idx[type][1]), 0));

  return results;
}
  
template <typename index_type, typename chirality_type>
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d<index_type, chirality_type>::
links_face_edge(index_type id)
{
  std::vector<std::pair<index_type, chirality_type> > results;
  
  int idx[4];
  fid2fidx(id, idx);
  // fprintf(stderr, "fid=%lu, fidx={%d, %d, %d, %d}, valid=%d\n",
  //     id, idx[0], idx[1], idx[2], idx[3], valid_fidx(idx));
  if (!valid_fidx(idx)) return results;
  
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
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d<index_type, chirality_type>::
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
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d<index_type, chirality_type>::
links_node_edge(index_type i)
{
  std::vector<std::pair<index_type, chirality_type> > results;
  // TODO
  return results;
}

template <typename index_type, typename chirality_type>
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d<index_type, chirality_type>::
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

  for (int p=0; p<4; p++)
    if (valid_fidx(contained_faces_fidx[t][p]))
      results.push_back(std::make_pair(fidx2fid(contained_faces_fidx[t][p]), contained_faces_chi[t][p]));

  return results;
}

template <typename index_type, typename chirality_type>
std::vector<std::pair<index_type, chirality_type> > mesh_graph_regular_3d<index_type, chirality_type>::links_face_cell(index_type id)
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
  for (int p=0; p<2; p++) 
    results.push_back(std::make_pair(cidx2cid(contained_cells_cidx[t][p]), contained_cells_chi[p]));

  return results;
}

template <typename index_type, typename chirality_type>
index_type mesh_graph_regular_3d<index_type, chirality_type>::nidx2nid(const int idx_[3]) const
{
  int idx[3] = {idx_[0], idx_[1], idx_[2]};
  for (int i=0; i<3; i++) {
    idx[i] = idx[i] % d[i];
    if (idx[i] < 0)
      idx[i] += d[i];
  }
  return idx[0] + d[0] * (idx[1] + d[1] * idx[2]); 
}

template <typename index_type, typename chirality_type>
index_type mesh_graph_regular_3d<index_type, chirality_type>::eidx2eid(const int idx[4]) const
{
  return nidx2nid(idx)*3 + idx[3];
}

template <typename index_type, typename chirality_type>
index_type mesh_graph_regular_3d<index_type, chirality_type>::fidx2fid(const int idx[4]) const
{
  return nidx2nid(idx)*3 + idx[3];
}

template <typename index_type, typename chirality_type>
index_type mesh_graph_regular_3d<index_type, chirality_type>::cidx2cid(const int idx[3]) const
{
  return nidx2nid(idx);
}

template <typename index_type, typename chirality_type>
bool mesh_graph_regular_3d<index_type, chirality_type>::valid_nidx(const int idx[3]) const
{
  for (int i=0; i<3; i++) 
    if (idx[i]<0 || idx[i]>=d[i]) 
      return false;
  return true;
}

template <typename index_type, typename chirality_type>
void mesh_graph_regular_3d<index_type, chirality_type>::nid2nidx(index_type id, int idx[3]) const
{
  int s = d[0] * d[1]; 
  int k = id / s; 
  int j = (id - k*s) / d[0]; 
  int i = id - k*s - j*d[0]; 

  idx[0] = i; idx[1] = j; idx[2] = k;
}

template <typename index_type, typename chirality_type>
void mesh_graph_regular_3d<index_type, chirality_type>::eid2eidx(index_type id, int idx[4]) const
{
  unsigned int nid = id / 3;
  nid2nidx(nid, idx);
  idx[3] = id % 3;
}

template <typename index_type, typename chirality_type>
void mesh_graph_regular_3d<index_type, chirality_type>::fid2fidx(index_type id, int idx[4]) const
{
  unsigned int nid = id / 3;
  nid2nidx(nid, idx);
  idx[3] = id % 3;
}

template <typename index_type, typename chirality_type>
void mesh_graph_regular_3d<index_type, chirality_type>::cid2cidx(index_type id, int idx[3]) const
{
  nid2nidx(id, idx);
}

template <typename index_type, typename chirality_type>
bool mesh_graph_regular_3d<index_type, chirality_type>::
valid_eidx(const int eidx[4]) const
{
  // fprintf(stderr, "eidx=%d, %d, %d, %d\n", eidx[0], eidx[1], eidx[2], eidx[3]);
  if (eidx[3]<0 || eidx[3]>=3) return false;
  else {
    int o[3] = {0};
    for (int i=0; i<3; i++) 
      if (pbc[i]) {
        if (eidx[i]<0 || eidx[i]>=d[i]) return false;
      } else {
        // if (eidx[i]<0 || eidx[i]>=d[i]-1) return false;
        if (eidx[i]<0 || eidx[i]>d[i]-1) return false;
        else if (eidx[i] == d[i]-1) o[i] = 1;
      }

    const int sum = o[0] + o[1] + o[2];
    if (sum == 0) return true;
    else if (sum > 2) return false;
    else if (o[0] && eidx[3] == 0) return false; 
    else if (o[1] && eidx[3] == 1) return false;
    else if (o[2] && eidx[3] == 2) return false;
    else return true;
    // return true;
  }
}

template <typename index_type, typename chirality_type>
bool mesh_graph_regular_3d<index_type, chirality_type>::
valid_fidx(const int fidx[4]) const
{
  if (fidx[3]<0 || fidx[3]>=3) return false;
  else {
    int o[3] = {0}; 
    for (int i=0; i<3; i++) 
      if (pbc[i]) {
        if (fidx[i]<0 || fidx[i]>=d[i]) return false;
      } else {
        if (fidx[i]<0 || fidx[i]>d[i]-1) return false;
        else if (fidx[i] == d[i]-1) o[i] = 1;
      }

    const int sum = o[0] + o[1] + o[2];
    if (sum == 0) return true;
    else if (o[0] + o[1] + o[2] > 1) return false;
    else if (o[0] && fidx[3] == 0) return true; 
    else if (o[1] && fidx[3] == 1) return true;
    else if (o[2] && fidx[3] == 2) return true;
    else return false;
  }
}

template <typename index_type, typename chirality_type>
bool mesh_graph_regular_3d<index_type, chirality_type>::valid_cidx(const int idx[3]) const
{
  for (int i=0; i<3; i++)
    if (pbc[i]) {
      if (idx[i] < 0 || idx[i] >= d[i]) return false;
    } else {
      if (idx[i] < 0 || idx[i] >= d[i]-1) return false;
    }
  return true;
}

}

#endif
