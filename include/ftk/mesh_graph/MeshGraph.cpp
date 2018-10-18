#include "MeshGraph.h"
#include <cassert>
#include <climits>
#include <cstdio>
#include <climits>

#if WITH_PROTOBUF
#include "MeshGraph.pb.h"
#include <google/protobuf/io/coded_stream.h> // for parsing large message
#endif

using std::make_tuple;
using std::get;

EdgeIdType2 AlternateEdge(EdgeIdType2 e, ChiralityType chirality)
{
  if (chirality>0)
    return e;
  else 
    return make_tuple(get<1>(e), get<0>(e));
}

FaceIdType3 AlternateFace(FaceIdType3 f, int rotation, ChiralityType chirality)
{
  if (chirality>0) {
    switch (rotation) {
    case 0: return f; 
    case 1: return make_tuple(get<2>(f), get<0>(f), get<1>(f));
    case 2: return make_tuple(get<1>(f), get<2>(f), get<0>(f));
    default: assert(false);
    }
  } else {
    switch (rotation) {
    case 0: return make_tuple(get<2>(f), get<1>(f), get<0>(f));
    case 1: return make_tuple(get<0>(f), get<2>(f), get<1>(f));
    case 2: return make_tuple(get<1>(f), get<0>(f), get<2>(f));
    default: assert(false);
    }
  }

  return make_tuple(UINT_MAX, UINT_MAX, UINT_MAX); // make compiler happy
}

FaceIdType4 AlternateFace(FaceIdType4 f, int rotation, ChiralityType chirality)
{
  using namespace std;

  if (chirality>0) {
    switch (rotation) {
    case 0: return f;
    case 1: return make_tuple(get<3>(f), get<0>(f), get<1>(f), get<2>(f));
    case 2: return make_tuple(get<2>(f), get<3>(f), get<0>(f), get<1>(f));
    case 3: return make_tuple(get<1>(f), get<2>(f), get<3>(f), get<0>(f));
    default: assert(false);
    }
  } else {
    switch (rotation) {
    case 0: return make_tuple(get<3>(f), get<2>(f), get<1>(f), get<0>(f));
    case 1: return make_tuple(get<0>(f), get<3>(f), get<2>(f), get<1>(f));
    case 2: return make_tuple(get<1>(f), get<0>(f), get<3>(f), get<2>(f));
    case 3: return make_tuple(get<2>(f), get<1>(f), get<0>(f), get<3>(f));
    default: assert(false);
    }
  }

  return make_tuple(UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX); // make compiler happy
}

////////////////////////
MeshGraph::~MeshGraph()
{
  Clear();
}

void MeshGraph::Clear()
{
  edges.clear();
  faces.clear();
  cells.clear();
}

void MeshGraph::SerializeToString(std::string &str) const
{
#if WITH_PROTOBUF
  PBMeshGraph pmg;

  for (int i=0; i<edges.size(); i++) {
    PBEdge *pedge = pmg.add_edges();
    pedge->set_node0( edges[i].node0 );
    pedge->set_node1( edges[i].node1 );

    for (int j=0; j<edges[i].contained_faces.size(); j++) {
      pedge->add_contained_faces( edges[i].contained_faces[j] );
      pedge->add_contained_faces_chirality( edges[i].contained_faces_chirality[j] );
      pedge->add_contained_faces_eid( edges[i].contained_faces_eid[j] );
    }
  }

  for (int i=0; i<faces.size(); i++) {
    PBFace *pface = pmg.add_faces();

    for (int j=0; j<faces[i].nodes.size(); j++) 
      pface->add_nodes(faces[i].nodes[j]);

    for (int j=0; j<faces[i].edges.size(); j++) {
      pface->add_edges(faces[i].edges[j]);
      pface->add_edges_chirality(faces[i].edges_chirality[j]);
    }

    for (int j=0; j<faces[i].contained_cells.size(); j++) {
      pface->add_contained_cells(faces[i].contained_cells[j]);
      pface->add_contained_cells_chirality(faces[i].contained_cells_chirality[j]);
      pface->add_contained_cells_fid(faces[i].contained_cells_fid[j]);
    }
  }

  for (int i=0; i<cells.size(); i++) {
    PBCell *pcell = pmg.add_cells();

    for (int j=0; j<cells[i].nodes.size(); j++)
      pcell->add_nodes(cells[i].nodes[j]);

    for (int j=0; j<cells[i].faces.size(); j++) {
      pcell->add_faces(cells[i].faces[j]);
      pcell->add_faces_chirality(cells[i].faces_chirality[j]);
      pcell->add_neighbor_cells(cells[i].neighbor_cells[j]);
    }
  }

  pmg.SerializeToString(&str);
#endif
}

bool MeshGraph::ParseFromString(const std::string& str)
{
#if WITH_PROTOBUF
  PBMeshGraph pmg;

  google::protobuf::io::CodedInputStream stream((uint8_t*)str.data(), str.size());
  stream.SetTotalBytesLimit(INT_MAX, INT_MAX/2);
  // if (!pmg.ParseFromCodedStream(&stream)) 
  if (!pmg.MergeFromCodedStream(&stream)) 
    return false;

  Clear();

  for (int i=0; i<pmg.edges_size(); i++) {
    PBEdge pedge = pmg.edges(i);
    CEdge edge;

    edge.node0 = pedge.node0();
    edge.node1 = pedge.node1();

    for (int j=0; j<pedge.contained_faces_size(); j++) {
      edge.contained_faces.push_back( pedge.contained_faces(j) );
      edge.contained_faces_chirality.push_back( pedge.contained_faces_chirality(j) );
      edge.contained_faces_eid.push_back( pedge.contained_faces_eid(j) );
    }

    edges.push_back(edge);
  }

  for (int i=0; i<pmg.faces_size(); i++) {
    PBFace pface = pmg.faces(i);
    CFace face;

    for (int j=0; j<pface.nodes_size(); j++)
      face.nodes.push_back( pface.nodes(j) );

    for (int j=0; j<pface.edges_size(); j++) {
      face.edges.push_back( pface.edges(j) );
      face.edges_chirality.push_back( pface.edges_chirality(j) );
    }

    for (int j=0; j<pface.contained_cells_size(); j++) {
      face.contained_cells.push_back( pface.contained_cells(j) );
      face.contained_cells_chirality.push_back( pface.contained_cells_chirality(j) );
      face.contained_cells_fid.push_back( pface.contained_cells_fid(j) );
    }

    faces.push_back(face);
  }

  for (int i=0; i<pmg.cells_size(); i++) {
    PBCell pcell = pmg.cells(i); 
    CCell cell;
  
    for (int j=0; j<pcell.nodes_size(); j++)
      cell.nodes.push_back( pcell.nodes(j) );

    for (int j=0; j<pcell.faces_size(); j++) {
      cell.faces.push_back( pcell.faces(j) );
      cell.faces_chirality.push_back( pcell.faces_chirality(j) );
    }

    for (int j=0; j<pcell.neighbor_cells_size(); j++)
      cell.neighbor_cells.push_back( pcell.neighbor_cells(j) );

    cells.push_back(cell);
  }
  return true;
#else
  return false;
#endif
}

void MeshGraph::SerializeToFile(const std::string& filename) const
{
  FILE *fp = fopen(filename.c_str(), "wb");

  std::string buf;
  SerializeToString(buf);
  fwrite(buf.data(), 1, buf.size(), fp);

  fclose(fp);
}

bool MeshGraph::ParseFromFile(const std::string& filename)
{
  FILE *fp = fopen(filename.c_str(), "rb"); 
  if (!fp) return false;

  fseek(fp, 0L, SEEK_END);
  size_t sz = ftell(fp);

  std::string buf;
  buf.resize(sz);
  fseek(fp, 0L, SEEK_SET);
  size_t sz1 = fread((char*)buf.data(), 1, sz, fp);
  fclose(fp);

  return ParseFromString(buf);
}

MeshGraphBuilder::MeshGraphBuilder(MeshGraph& mg)
  : _mg(mg)
{
}


////////////////////////
////////////////
MeshGraphBuilder_Tet::MeshGraphBuilder_Tet(int ncells, MeshGraph& mg) : 
  MeshGraphBuilder(mg) 
{
  mg.cells.resize(ncells);
}

EdgeIdType MeshGraphBuilder_Tet::GetEdge(EdgeIdType2 e2, ChiralityType &chirality)
{
  for (chirality=-1; chirality<2; chirality+=2) {
    std::map<EdgeIdType2, EdgeIdType>::iterator it = _edge_map.find(AlternateEdge(e2, chirality)); 
    if (it != _edge_map.end())
      return it->second;
  }
  return UINT_MAX;
}

FaceIdType MeshGraphBuilder_Tet::GetFace(FaceIdType3 f3, ChiralityType &chirality)
{
  for (chirality=-1; chirality<2; chirality+=2) 
    for (int rotation=0; rotation<3; rotation++) {
      std::map<FaceIdType3, FaceIdType>::iterator it = _face_map.find(AlternateFace(f3, rotation, chirality));
      if (it != _face_map.end())
        return it->second;
    }
  return UINT_MAX;
}

EdgeIdType MeshGraphBuilder_Tet::AddEdge(EdgeIdType2 e2, ChiralityType &chirality, FaceIdType f, int eid)
{
  EdgeIdType e = GetEdge(e2, chirality);

  if (e == UINT_MAX) {
    e = _mg.edges.size();
    _edge_map.insert(std::pair<EdgeIdType2, EdgeIdType>(e2, e));
    
    CEdge edge1;
    edge1.node0 = get<0>(e2);
    edge1.node1 = get<1>(e2);
    
    _mg.edges.push_back(edge1);
    chirality = 1;
  }
  
  CEdge &edge = _mg.edges[e];

  edge.contained_faces.push_back(f);
  edge.contained_faces_chirality.push_back(chirality);
  edge.contained_faces_eid.push_back(eid);
  
  return e;
}

FaceIdType MeshGraphBuilder_Tet::AddFace(FaceIdType3 f3, ChiralityType &chirality, CellIdType c, int fid)
{
  FaceIdType f = GetFace(f3, chirality);

  if (f == UINT_MAX) {
    f = _face_map.size();
    _face_map.insert(std::pair<FaceIdType3, FaceIdType>(f3, f));
    
    CFace face1;
    face1.nodes.push_back(get<0>(f3));
    face1.nodes.push_back(get<1>(f3));
    face1.nodes.push_back(get<2>(f3));

    EdgeIdType2 e2[3] = {
      make_tuple(get<0>(f3), get<1>(f3)),
      make_tuple(get<1>(f3), get<2>(f3)),
      make_tuple(get<2>(f3), get<0>(f3))};

    for (int i=0; i<3; i++) {
      EdgeIdType e = AddEdge(e2[i], chirality, f, i);
      face1.edges.push_back(e);
      face1.edges_chirality.push_back(chirality);
    }
    
    _mg.faces.push_back(face1);
    chirality = 1;
  }
  
  CFace &face = _mg.faces[f];

  face.contained_cells.push_back(c);
  face.contained_cells_chirality.push_back(chirality);
  face.contained_cells_fid.push_back(fid);

  return f;
}

void MeshGraphBuilder_Tet::AddCell(
    CellIdType c,
    const std::vector<NodeIdType> &nodes, 
    const std::vector<CellIdType> &neighbors, 
    const std::vector<FaceIdType3> &faces)
{
  CCell &cell = _mg.cells[c];

  // nodes
  cell.nodes = nodes;

  // neighbor cells
  cell.neighbor_cells = neighbors;

  // faces and edges
  for (int i=0; i<faces.size(); i++) {
    ChiralityType chirality; 
    FaceIdType3 f3 = faces[i];

    FaceIdType f = AddFace(f3, chirality, c, i);
    cell.faces.push_back(f);
    cell.faces_chirality.push_back(chirality);
  }
}
