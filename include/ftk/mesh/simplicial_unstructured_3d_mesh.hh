#ifndef _FTK_MESH_SIMPLEX_3D_HH
#define _FTK_MESH_SIMPLEX_3D_HH

#include <ftk/ftk_config.hh>
#include <ftk/object.hh>
#include <ftk/ndarray.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/algorithms/bfs.hh>
#include <ftk/external/diy-ext/serialization.hh>
#include <set>
#include <iostream>
#include <vector>

#if FTK_HAVE_VTK
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnsignedLongArray.h>
#include <vtkLongArray.h>
#include <vtkGenericCell.h>
#include <vtkDataSetWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkPointData.h>
#include <vtkPoints2D.h>
#endif

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_unstructured_3d_mesh : public object {
  simplicial_unstructured_3d_mesh() {}

  simplicial_unstructured_3d_mesh(
      const std::vector<F>& coords, 
      const std::vector<I>& tetrahedra);

  size_t nd() const {return 3;}

  void build_edges();
  void build_triangles();
  void build_tetrahedra();

public: // io
  void from_vtk_unstructured_grid_file(const std::string &filename);
  void to_vtk_unstructured_grid_file(const std::string& filename) const;

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
inline void simplicial_unstructured_3d_mesh<I, F>::from_vtk_unstructured_grid_file(const std::string &filename)
{
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkXMLUnstructuredGridReader::New();
  reader->SetFileName(filename.c_str());
  reader->Update();
  from_vtk_unstructured_grid(reader->GetOutput());
#else
  fatal("FTK not commpiled with VTK");
#endif
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
#endif

}

#endif
