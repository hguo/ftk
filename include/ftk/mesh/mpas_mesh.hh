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

  void initialize_c2v_interpolants();
  ndarray<F> interpolate_velocity_c2v(const ndarray<F>& Vc) const;

#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> surface_cells_to_vtu() const;
#endif
  void surface_cells_to_vtu(const std::string filename) const;

public:
  size_t n_cells() const { return xyzCells.dim(1); }
  size_t n_vertices() const { return xyzVertices.dim(1); }

  size_t locate_cell_i(const F x[]) const { return kd_cells->find_nearest(x); }
  bool point_in_cell(const F x[]) const; // WIP
  
  // void interpolate_c2v(const I vid, const F cvals[3][], F vvals[]) const;

  I cid2i(I cid) const { auto it = cellIdToIndex.find(cid); if (it != cellIdToIndex.end()) return it->second; else return -1;  }
  I vid2i(I vid) const { auto it = vertexIdToIndex.find(vid); if (it != vertexIdToIndex.end()) return it->second; else return -1; }

  I i2cid(I i) const { return indexToCellID[i]; }
  I i2vid(I i) const { return indexToVertexID[i]; }

  I verts_i_on_cell_i(const I i, I vi[]) const; // returns number of verts
  void vert_i_coords(const I i, F X[3]) const;
  void verts_i_coords(const I n, const I i[], F X[][3]) const;

  bool is_vertex_on_boundary(I vid) const { return vertex_on_boundary[vertex_id_to_index(vid)]; }

public:
  std::shared_ptr<kd_t<F, 3>> kd_cells, kd_vertices;
  
  ndarray<I> cellsOnVertex, verticesOnCell;
  ndarray<I> indexToVertexID, indexToEdgeID, indexToCellID;
  ndarray<F> xyzCells, xyzVertices;
  ndarray<I> nEdgesOnCell;

  ndarray<F> c2v_interpolants;
  std::vector<bool> vertex_on_boundary;

  std::map<I, I> cellIdToIndex, // cellID->index
                 edgeIdToIndex,              
                 vertexIdToIndex; // vertexID->index

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
  kd_cells.reset(new kd_t<F, 3>);
  kd_cells->set_inputs(this->xyzCells);
  kd_cells->build();
}

template <typename I, typename F>
void mpas_mesh<I, F>::initialize_c2v_interpolants()
{
  c2v_interpolants.reshape(3, n_vertices());
  vertex_on_boundary.resize(n_vertices());

  for (I i = 0; i < n_vertices(); i ++) {
    F Xc[3][3], Xv[3];
   
    bool boundary = false;
    for (I j = 0; j < 3; j ++) {
      const I c = cid2i( cellsOnVertex(j, i) );
      if (c < 0) {
        boundary = true;
        break;
      }
      // fprintf(stderr, "cell=%d\n", c);
      Xv[j] = xyzVertices(j, i);
      for (I k = 0; k < 3; k ++)
        Xc[j][k] = xyzCells(k, c);
    }

    if (boundary) {
      vertex_on_boundary[i] = true;
      continue;
    }

    // print3x3("Xc", Xc);
    // print3("Xv", Xv);
      
    F mu[3];
    bool succ = inverse_lerp_s2v3_0(Xc, Xv, mu);
    for (I k = 0; k < 3; k ++)
      c2v_interpolants(k, i) = mu[k];

    // fprintf(stderr, "mu=%f, %f, %f\n", 
    //     mu[0], mu[1], mu[2]);
  }
}

template <typename I, typename F>
ndarray<F> mpas_mesh<I, F>::interpolate_velocity_c2v(const ndarray<F>& Vc) const
{
  ndarray<F> V; // velocities on vertices
  V.reshape({3, Vc.dim(1), n_vertices()});

  for (auto i = 0; i < n_vertices(); i ++) {
    const auto vid = i2vid(i);
    // bool boundary = is_vertex_on_boundary(vid);
    const bool boundary = vertex_on_boundary[i];

    F lambda[3] = {0};
    I c[3] = {0};
    if (!boundary) {
      for (auto k = 0; k < 3; k ++) {
        lambda[k] = c2v_interpolants(k, i);
        c[k] = cid2i( cellsOnVertex(k, i) );
      }
    }

    for (auto layer = 0; layer < Vc.dim(1); layer ++) {
      for (auto k = 0; k < 3; k ++) {
        if (boundary) {
          V(k, layer, i) = F(0);
        } else {
          for (auto l = 0; l < 3; l ++) 
            V(k, layer, i) += lambda[l] * Vc(k, layer, c[l]);
        }
      }
    }
  }

  return V;
}

// inline mpas_mesh::mpas_mesh(const mpas_mesh& mm) :

template <typename I, typename F>
void mpas_mesh<I, F>::vert_i_coords(const I i, F X[3]) const
{
  for (int k = 0; k < 3; k ++)
    X[k] = xyzVertices(k, i);
}

template <typename I, typename F>
void mpas_mesh<I, F>::verts_i_coords(const I n, const I i[], F X[][3]) const
{
  for (auto j = 0; j < n; j ++)
    vert_i_coords(i[j], X[j]);
}

template <typename I, typename F>
I mpas_mesh<I, F>::verts_i_on_cell_i(const I i, I vi[]) const
{
  I n = nEdgesOnCell[i];
  for (auto j = 0; j < n; j ++)
    vi[j] = vid2i( verticesOnCell(j, i) );
  return n;
}

template <typename I, typename F>
void mpas_mesh<I, F>::read_netcdf(const std::string filename, diy::mpi::communicator comm)
{
#if FTK_HAVE_NETCDF
  int ncid;
  NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );

  cellsOnVertex.read_netcdf(ncid, "cellsOnVertex"); // "conn"
  verticesOnCell.read_netcdf(ncid, "verticesOnCell");
  nEdgesOnCell.read_netcdf(ncid, "nEdgesOnCell");

  indexToCellID.read_netcdf(ncid, "indexToCellID");
  for (auto i = 0; i < indexToCellID.size(); i ++)
    cellIdToIndex[indexToCellID[i]] = i;

  indexToEdgeID.read_netcdf(ncid, "indexToEdgeID");
  for (auto i = 0; i < indexToEdgeID.size(); i ++)
    edgeIdToIndex[indexToEdgeID[i]] = i;
  
  indexToVertexID.read_netcdf(ncid, "indexToVertexID");
  for (auto i = 0; i < indexToVertexID.size(); i ++) 
    vertexIdToIndex[indexToVertexID[i]] = i;

  std::vector<I> vconn;
  for (int i = 0; i < cellsOnVertex.dim(1); i ++) {
    if (cellsOnVertex(0, i) == 0 || 
        cellsOnVertex(1, i) == 0 ||
        cellsOnVertex(2, i) == 0) 
      continue;
    else {
      for (int j = 0; j < 3; j ++) 
        vconn.push_back(cellIdToIndex[cellsOnVertex(j, i)]);
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
#endif
}

#if FTK_HAVE_VTK
template <typename I, typename F>
vtkSmartPointer<vtkUnstructuredGrid> mpas_mesh<I, F>::surface_cells_to_vtu() const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> pts = vtkPoints::New();
  pts->SetNumberOfPoints(n_vertices());

  for (auto i = 0; i < n_vertices(); i ++)
    pts->SetPoint(i, xyzVertices(0, i), xyzVertices(1, i), xyzVertices(2, i));
  grid->SetPoints(pts);

  for (auto i = 0; i < n_cells(); i ++) {
    vtkIdType ids[7];
    for (auto j = 0; j < 7; j ++) {
      auto id = vid2i( verticesOnCell(j, i) );
      ids[j] = id;
    }

    const int n = nEdgesOnCell[i];
    // fprintf(stderr, "%zu, %d\n", i, n);
      
    grid->InsertNextCell( VTK_POLYGON, n, ids );

#if 0
    if (n == 6)
      grid->InsertNextCell( VTK_HEXAHEDRON, n, ids );
    else if (n == 7 || n == 5) 
      grid->InsertNextCell( VTK_POLYGON, n, ids );
    else if (n == 4)
      grid->InsertNextCell( VTK_QUAD, n, ids );
#endif
  }

  return grid;
}
#endif
  
template <typename I, typename F>
void mpas_mesh<I, F>::surface_cells_to_vtu(const std::string filename) const
{
#if FTK_HAVE_VTK
  auto grid = surface_cells_to_vtu();
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInputData( grid );
  writer->Write();
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

} // namespace ftk

#endif 
