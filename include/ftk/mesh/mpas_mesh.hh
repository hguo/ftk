#ifndef _FTK_MPAS_2D_MESH_HH
#define _FTK_MPAS_2D_MESH_HH

#include <ftk/config.hh>
#include <ftk/basic/kd.hh>
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct mpas_mesh { // : public simplicial_unstructured_2d_mesh<I, F> {
  void read_netcdf(const std::string filename, diy::mpi::communicator comm = MPI_COMM_WORLD);
  void initialize();

public:
  size_t nCells() const { return xyzCells.dim(1); }
  size_t nVertices() const { return xyzVertices.dim(1); }

  size_t locateCell(const std::array<F, 3> x) const { return kdCells->find_nearest(x); }

public:
  std::shared_ptr<kd_t<F, 3>> kdCells, kdVertices;
  
  ndarray<I> cellsOnVertex, verticesOnCell;
  ndarray<I> indexToVertexID, indexToEdgeID, indexToCellID;
  ndarray<F> xyzCells, xyzVertices;

#if 0
  ndarray<I> voronoi_c2v, // (6, nc)
             voronoi_v2c; // (3, nv)
  ndarray<F> voronoi_vcoords; // (3, nv)
#endif
};

//////
template <typename I, typename F>
void mpas_mesh<I, F>::initialize()
{
  kdCells.reset(new kd_t<F, 3>);
  kdCells->set_inputs(this->xyzVertices);
  kdCells->build();
}

// inline mpas_mesh::mpas_mesh(const mpas_mesh& mm) :

template <typename I, typename F>
void mpas_mesh<I, F>::read_netcdf(const std::string filename, diy::mpi::communicator comm)
{
#if FTK_HAVE_NETCDF
  int ncid;
  NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );

  cellsOnVertex.read_netcdf(ncid, "cellsOnVertex"); // "conn"
  verticesOnCell.read_netcdf(ncid, "verticesOnCell");
  indexToCellID.read_netcdf(ncid, "indexToCellID");

  std::map<I, I> cellIndex;
  for (int i = 0; i < indexToCellID.size(); i ++)
    cellIndex[indexToCellID[i]] = i;

  std::vector<I> vconn;
  for (int i = 0; i < cellsOnVertex.dim(1); i ++) {
    if (cellsOnVertex(0, i) == 0 || 
        cellsOnVertex(1, i) == 0 ||
        cellsOnVertex(2, i) == 0) 
      continue;
    else {
      for (int j = 0; j < 3; j ++) 
        vconn.push_back(cellIndex[cellsOnVertex(j, i)]);
    }
  }

  ndarray<I> conn;
  conn.copy_vector(vconn);
  conn.reshape(3, vconn.size()/3);
 
  ndarray<F> xCell, yCell, zCell;
  xCell.read_netcdf(ncid, "xCell");
  yCell.read_netcdf(ncid, "yCell");
  zCell.read_netcdf(ncid, "zCell");
  // latCell.from_netcdf(ncid, "latCell");
  // lonCell.from_netcdf(ncid, "lonCell");
  
  xyzCells.reshape(3, xCell.size());
  for (size_t i = 0; i < xCell.size(); i ++) {
    xyzCells(0, i) = xCell[i];
    xyzCells(1, i) = yCell[i];
    xyzCells(2, i) = zCell[i];
  }

  ndarray<F> xVertex, yVertex, zVertex;
  xVertex.read_netcdf(ncid, "xVertex");
  yVertex.read_netcdf(ncid, "yVertex");
  zVertex.read_netcdf(ncid, "zVertex");

  xyzVertices.reshape(3, xVertex.size());
  for (size_t i = 0; i < xVertex.size(); i ++) {
    xyzVertices(0, i) = xVertex[i];
    xyzVertices(1, i) = yVertex[i];
    xyzVertices(2, i) = zVertex[i];
  }

  NC_SAFE_CALL( nc_close(ncid) );

  fprintf(stderr, "mpas mesh: #vert=%zu, #tri=%zu\n", 
      xyzCells.dim(1), conn.dim(1));

  // return std::shared_ptr<mpas_mesh<I, F>>(
  //     new mpas_mesh<I, F>(xyz, conn));
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
  return NULL;
#endif
}

} // namespace ftk

#endif 
