#ifndef _FTK_MPAS_2D_MESH_HH
#define _FTK_MPAS_2D_MESH_HH

#include <ftk/config.hh>
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_mpas_2d_mesh : public simplicial_unstructured_2d_mesh<I, F> {
  simplicial_mpas_2d_mesh(
      const ndarray<F>& coords, 
      const ndarray<I>& triangles) : simplicial_unstructured_2d_mesh<I, F>(coords, triangles) {}

  // int ncoords() const { return 3; }

  static std::shared_ptr<simplicial_mpas_2d_mesh<I, F>> from_file(const std::string& filename, diy::mpi::communicator comm = MPI_COMM_WORLD);
};

// inline simplicial_mpas_2d_mesh::simplicial_mpas_2d_mesh(const mpas_mesh& mm) :

template <typename I, typename F>
std::shared_ptr<simplicial_mpas_2d_mesh<I, F>> simplicial_mpas_2d_mesh<I, F>::from_file(
    const std::string& filename, 
    diy::mpi::communicator comm)
{
#if FTK_HAVE_NETCDF
  int ncid;
  NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );

  ndarray<I> cellsOnVertex, indexToCellID;
  cellsOnVertex.read_netcdf(ncid, "cellsOnVertex"); // "conn"
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
  
  ndarray<F> xCell, yCell, zCell, xyz;
  xCell.read_netcdf(ncid, "xCell");
  yCell.read_netcdf(ncid, "yCell");
  zCell.read_netcdf(ncid, "zCell");
  // latCell.from_netcdf(ncid, "latCell");
  // lonCell.from_netcdf(ncid, "lonCell");
  xyz.reshape(3, xCell.size());
  
  for (int i = 0; i < xCell.size(); i ++) {
    xyz(0, i) = xCell[i];
    xyz(1, i) = yCell[i];
    xyz(2, i) = zCell[i];
  }

  NC_SAFE_CALL( nc_close(ncid) );

  fprintf(stderr, "mpas mesh: #vert=%zu, #tri=%zu\n", 
      xyz.dim(1), conn.dim(1));

  return std::shared_ptr<simplicial_mpas_2d_mesh<I, F>>(
      new simplicial_mpas_2d_mesh<I, F>(xyz, conn));
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
  return NULL;
#endif
}

} // namespace ftk

#endif 
