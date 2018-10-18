#ifndef _FTK_MESHGRAPH_H
#define _FTK_MESHGRAPH_H

#include <vector>
#include <bitset>
#include <map>

#include <tuple>

namespace ftk {

struct CEdge;
struct CFace;
struct CCell;

typedef std::tuple<NodeIdType, NodeIdType> EdgeIdType2;
typedef std::tuple<NodeIdType, NodeIdType, NodeIdType> FaceIdType3;
typedef std::tuple<NodeIdType, NodeIdType, NodeIdType, NodeIdType> FaceIdType4;

EdgeIdType2 AlternateEdge(EdgeIdType2 e2, ChiralityType chirality);
FaceIdType3 AlternateFace(FaceIdType3 f3, int rotation, ChiralityType chirality);
FaceIdType4 AlternateFace(FaceIdType4 f4, int rotation, ChiralityType chirality);

struct CEdge {
  // nodes
  NodeIdType node0, node1;

  // neighbor faces (unordered)
  std::vector<FaceIdType> contained_faces;
  std::vector<ChiralityType> contained_faces_chirality;
  std::vector<int> contained_faces_eid; // the edge id in the corresponding face

  CEdge() : node0(0), node1(0) {
    // contained_faces.reserve(12);
    // contained_faces_chirality.reserve(12);
    // contained_faces_eid.reserve(12);
  }

  bool Valid() const {return node0 != node1;}
};

struct CFace {
  // nodes (ordered)
  std::vector<NodeIdType> nodes;

  // edges (ordered)
  std::vector<EdgeIdType> edges;
  std::vector<ChiralityType> edges_chirality;

  // neighbor cells, only two, chirality(cell0)=-1, chirality(cell1)=1
  // CellIdType contained_cell0, contained_cell1;
  // int contained_cell0_fid, contained_cell1_fid;
  std::vector<FaceIdType> contained_cells;
  std::vector<ChiralityType> contained_cells_chirality;
  std::vector<int> contained_cells_fid;

  // utils
  // bool on_boundary() const {return contained_cell0 == NULL || contained_cell1 == NULL;}
  CFace() {
    // nodes.reserve(3);
    // edges.reserve(6);
    // edges_chirality.reserve(6);
    // contained_cells.reserve(2);
    // contained_cells_chirality.reserve(2);
  }

  bool Valid() const {return nodes.size()>0;}
};

struct CCell {
  // nodes (ordered)
  std::vector<NodeIdType> nodes;

  // faces (ordered)
  std::vector<FaceIdType> faces;
  std::vector<ChiralityType> faces_chirality;

  // neighbor cells (ordered)
  std::vector<CellIdType> neighbor_cells;

  CCell() {
    // nodes.reserve(4);
    // faces.reserve(4);
    // faces_chirality.reserve(4);
    // neighbor_cells.reserve(4);
  }

  bool Valid() const {return faces.size()>0;}
};

class MeshGraphBuilder;
class MeshGraphBuilder_Tet;
class MeshGraphBuilder_Hex;

class MeshGraphRegular2D;
class MeshGraphRegular3D;

class MeshGraph {
protected:
  friend class MeshGraphBuilder;
  friend class MeshGraphBuilder_Tet;
  friend class MeshGraphBuilder_Hex;
  
  std::vector<CEdge> edges;
  std::vector<CFace> faces;
  std::vector<CCell> cells;

public:
  ~MeshGraph();
  
  void Clear();

  virtual EdgeIdType NEdges() const {return edges.size();}
  virtual FaceIdType NFaces() const {return faces.size();}
  virtual CellIdType NCells() const {return cells.size();}

  virtual CEdge Edge(EdgeIdType i, bool nodes_only=false) const {return edges[i];}
  virtual CFace Face(FaceIdType i, bool nodes_only=false) const {return faces[i];} // second arg for acceleration
  virtual CCell Cell(CellIdType i, bool nodes_only=false) const {return cells[i];}

  void SerializeToString(std::string &str) const;
  bool ParseFromString(const std::string &str);

  void SerializeToFile(const std::string& filename) const;
  bool ParseFromFile(const std::string& filename);
};

class MeshGraphBuilder {
public:
  explicit MeshGraphBuilder(MeshGraph& mg);
  virtual ~MeshGraphBuilder() {}
 
  // virtual void Build() = 0;

protected:
  MeshGraph &_mg;
};

class MeshGraphBuilder_Tet : public MeshGraphBuilder {
public:
  explicit MeshGraphBuilder_Tet(int ncells, MeshGraph& mg);
  ~MeshGraphBuilder_Tet() {}

  void AddCell(
      CellIdType c, 
      const std::vector<NodeIdType> &nodes, 
      const std::vector<CellIdType> &neighbors, 
      const std::vector<FaceIdType3> &faces);

  // void Build();

private:
  EdgeIdType AddEdge(EdgeIdType2 e2, ChiralityType &chirality, FaceIdType f, int eid);
  FaceIdType AddFace(FaceIdType3 f3, ChiralityType &chirality, CellIdType c, int fid);
  EdgeIdType GetEdge(EdgeIdType2 e2, ChiralityType &chirality);
  FaceIdType GetFace(FaceIdType3 f3, ChiralityType &chirality);

private:
  std::map<EdgeIdType2, EdgeIdType> _edge_map;
  std::map<FaceIdType3, FaceIdType> _face_map;
};

}

#endif
