#ifndef _FTK_MESH_SIMPLEX_3D_HH
#define _FTK_MESH_SIMPLEX_3D_HH

#include <ftk/ftk_config.hh>
#include <ftk/mesh/simplicial_unstructured_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_unstructured_3d_mesh : public simplicial_unstructured_mesh<I, F> {
  simplicial_unstructured_3d_mesh() {}

  simplicial_unstructured_3d_mesh(
      const std::vector<F>& coords, 
      const std::vector<I>& tetrahedra);

  int nd() const {return 3;}
  size_t n(int d) const;

  void build_edges();
  void build_triangles();
  void build_tetrahedra();

public: // io
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_vtk_unstructured_grid() const;
  void from_vtk_unstructured_grid(vtkSmartPointer<vtkUnstructuredGrid> grid);
#endif

public: 
  void element_for(int d, std::function<void(I)> f);

public:
  std::set<I> sides(int d, I i) const;
  std::set<I> side_of(int d, I i) const;

  void get_simplex(int d, I i, I simplex[]) const;
  void get_coords(I i, F coords[]) const;

  const ndarray<F>& get_coords() const {return vertex_coords;}

private: // mesh conn
  ndarray<F> vertex_coords; // 3*n_vert
  std::vector<std::set<I>> vertex_side_of;
  std::vector<std::set<I>> vertex_edge_vertex;

  ndarray<I> edges;
  ndarray<I> edges_side_of;

  ndarray<I> triangles;
  ndarray<I> triangle_sides;

  ndarray<I> tetrahedra;
  ndarray<I> tetrahedra_sides;

public:
  std::map<std::tuple<I, I>, int> edge_id_map;
  std::map<std::tuple<I, I, I>, int> triangle_id_map;
  std::map<std::tuple<I, I, I, I>, int> tetrahedron_id_map;
};

//////////
template <typename I, typename F>
size_t simplicial_unstructured_3d_mesh<I, F>::n(int d) const
{
  if (d == 0) 
    return vertex_coords.dim(1);
  else if (d == 1)
    return edges.dim(1);
  else if (d == 2)
    return triangles.dim(1);
  else if (d == 3)
    return tetrahedra.dim(1);
  else return 0;
}

#if FTK_HAVE_VTK  
template <typename I, typename F>
void simplicial_unstructured_3d_mesh<I, F>::from_vtk_unstructured_grid(vtkSmartPointer<vtkUnstructuredGrid> grid)
{
  vtkIdType ncells = grid->GetNumberOfCells();
  std::vector<int> tets;
  for (vtkIdType i = 0; i < ncells; i ++) {
    vtkSmartPointer<vtkCell> cell = grid->GetCell(i);
    if (cell->GetCellType() == VTK_TETRA) {
      vtkIdType v[4] = {
        cell->GetPointId(0), 
        cell->GetPointId(1), 
        cell->GetPointId(2),
        cell->GetPointId(3)
      };
      std::sort(v, v+4);
      for (int j = 0; j < 4; j ++)
        tets.push_back(v[j]);
    }
  }
  tetrahedra.reshape({4, tets.size()/4});
  tetrahedra.from_vector(tets);

  vtkIdType npts = grid->GetNumberOfPoints();
  vertex_coords.reshape({3, size_t(npts)});
  for (vtkIdType i = 0; i < npts; i ++) {
    double x[3];
    grid->GetPoint(i, x);
    for (int j = 0; j < 3; j ++)
      vertex_coords(j, i) = x[j];
  }

  // fprintf(stderr, "#tet=%zu, #pts=%zu\n", 
  //     tetrahedra.dim(1), vertex_coords.dim(1));

  // build_triangles();
  // build_edges();
}

template <typename I, typename F>
vtkSmartPointer<vtkUnstructuredGrid> simplicial_unstructured_3d_mesh<I, F>::
to_vtk_unstructured_grid() const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> pts = vtkPoints::New();
  pts->SetNumberOfPoints(n(0));

  for (int i=0; i<n(0); i++)
    pts->SetPoint(i, vertex_coords[i*3], vertex_coords[i*3+1], vertex_coords[i*3+2]); 

  for (int i=0; i<n(3); i ++) {
    vtkIdType ids[4] = {tetrahedra[i*4], tetrahedra[i*4+1], tetrahedra[i*4+2], tetrahedra[i*4+3]};
    grid->InsertNextCell(VTK_TETRA, 4, ids);
  }

  grid->SetPoints(pts);
  return grid;
}
#endif

}

#endif
