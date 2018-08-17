#ifndef _TRIANGULAR_MESH_2D
#define _TRIANGULAR_MESH_2D

#include "ftk/mesh/mesh.h"

namespace ftk {

struct TriangularMesh2D : public Mesh
{
  std::set<size_t> vertexNeighbors(size_t i) {return vertexNeighborTable[i];}

  size_t nVertices, nTriangles;
  std::vector<size_t> conn;
  // std::vector<double> coords;

  std::vector<std::set<size_t> > vertexNeighborTable;
};

}
