#include <ftk/mesh_graph/mesh_graph_regular_3d.hh>
#include <ftk/mesh_graph/mesh_graph_regular_3d_tets.hh>
#include <iostream>

int main(int argc, char **argv)
{
  ftk::mesh_graph_regular_3d_tets<> mg(128, 128, 128);
#if 0
  size_t eid = 889;
  int eidx[4];
  mg.eid2eidx(eid, eidx);
  const auto valid = mg.valid(1, eid);

  fprintf(stderr, "valid=%d, eidx={%d, %d, %d, %d}\n", 
      valid, eidx[0], eidx[1], eidx[2], eidx[3]);

  const auto nodes = mg.links_edge_node(eid);
  fprintf(stderr, "%lu\n", nodes.size());
#endif

  size_t fid = 902465; 
  int fidx[4];
  mg.fid2fidx(fid, fidx);
  fprintf(stderr, "fidx={%d, %d, %d, %d}\n", fidx[0], fidx[1], fidx[2], fidx[3]);
#if 0
  const auto cells = mg.links_face_cell(fid);
  fprintf(stderr, "#cells=%lu\n", cells.size());

  for (auto c : cells) {
    const auto faces = mg.links_cell_face(c.first);
    int cidx[4];
    mg.cid2cidx(c.first, cidx);
    fprintf(stderr, "-cid=%lu, cidx={%lu, %lu, %lu, %lu}, #faces=%lu\n", 
        c.first, cidx[0], cidx[1], cidx[2], cidx[3], faces.size());
    for (auto f : faces) {
      fprintf(stderr, "--fid=%lu\n", f.first);
    }
  }
#endif

  const auto edges = mg.links_face_edge(fid);
  for (auto e : edges) {
    int eidx[4];
    mg.eid2eidx(e.first, eidx);
    fprintf(stderr, "-eid=%lu, eidx={%d, %d, %d, %d}, chi=%d\n", e.first, eidx[0], eidx[1], eidx[2], eidx[3], e.second);
    const auto nodes = mg.links_edge_node(e.first);
    for (auto n : nodes) {
      int nidx[3];
      mg.nid2nidx(n.first, nidx);
      fprintf(stderr, "--nid=%lu, nidx={%d, %d, %d}\n", n.first, nidx[0], nidx[1], nidx[2]);
    }
  }

  fprintf(stderr, "----\n");
  const auto nodes = mg.links_face_node(fid); 
  for (auto n : nodes) {
    fprintf(stderr, "nid=%lu\n", n.first);
  }

#if 0
  const auto valid = mg.valid(2, fid);
  const auto faces = mg.links_face_edge(fid);
  const auto nodes = mg.links_face_node(fid);
  fprintf(stderr, "valid=%d, nfaces=%lu, nnodes=%lu\n", valid, faces.size(), nodes.size());
#endif

  return 0;
}
