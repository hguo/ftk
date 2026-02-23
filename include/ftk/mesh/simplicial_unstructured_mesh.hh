#ifndef _FTK_UNSTRUCTURED_MESH_HH
#define _FTK_UNSTRUCTURED_MESH_HH

#include <ftk/config.hh>
#include <ftk/object.hh>
#include <ndarray/ndarray.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/algorithms/bfs.hh>
#include <ftk/utils/serialization.hh>
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
#include <vtkCellTypes.h>
#endif

namespace ftk {

/**
 * @brief Base class for simplicial unstructured meshes
 *
 * This class provides the interface for unstructured meshes composed of
 * simplices (triangles in 2D, tetrahedra in 3D). It supports I/O with VTK
 * formats and provides methods for associating scalar and vector data with
 * mesh vertices.
 *
 * Simplicial meshes are commonly used in finite element analysis and
 * computational geometry. FTK uses them as the foundation for feature
 * tracking on unstructured grids.
 *
 * @tparam I Index type (typically int or long)
 * @tparam F Floating-point type for coordinates (typically float or double)
 */
template <typename I=int, typename F=double>
struct simplicial_unstructured_mesh : public object {
  simplicial_unstructured_mesh() {}

  /**
   * @brief Get the dimensionality of the mesh
   * @return Mesh dimension (2 for 2D triangular meshes, 3 for 3D tetrahedral meshes)
   */
  virtual int nd() const = 0;

  /**
   * @brief Get the number of d-dimensional elements
   * @param d Dimension of elements (0=vertices, 1=edges, 2=faces, 3=cells)
   * @param part If true, return count for local partition only
   * @return Number of d-dimensional elements
   */
  virtual size_t n(int d, bool part = false) const = 0;

public: // io
  /**
   * @brief Load mesh from legacy VTK file format
   * @param filename Path to .vtk file
   */
  void from_legacy_vtk_file(const std::string& filename);

  /**
   * @brief Load mesh from VTK XML unstructured grid file (.vtu)
   * @param filename Path to .vtu file
   */
  void from_vtk_unstructured_grid_file(const std::string &filename);

  /**
   * @brief Save mesh to VTK XML unstructured grid file (.vtu)
   * @param filename Path to output .vtu file
   */
  void to_vtk_unstructured_grid_file(const std::string &filename) const;

  /**
   * @brief Save mesh with scalar data to VTK file
   * @param filename Path to output .vtu file
   * @param varname Name of the scalar variable
   * @param scalar Scalar field data (one value per vertex)
   */
  void scalar_to_vtk_unstructured_grid_data_file(const std::string& filename, const std::string& varname, const ndarray<F>&) const;

  /**
   * @brief Save mesh with vector data to VTK file
   * @param filename Path to output .vtu file
   * @param varname Name of the vector variable
   * @param vector Vector field data (2 or 3 components per vertex)
   */
  void vector_to_vtk_unstructured_grid_data_file(const std::string& filename, const std::string& varname, const ndarray<F>&) const;
#if FTK_HAVE_VTK
  /**
   * @brief Create VTK unstructured grid with scalar data
   * @param varname Name of the scalar variable
   * @param scalar Scalar field data
   * @return VTK unstructured grid object
   */
  vtkSmartPointer<vtkUnstructuredGrid> scalar_to_vtk_unstructured_grid_data(const std::string& varname, const ndarray<F>&) const;

  /**
   * @brief Create VTK unstructured grid with vector data
   * @param varname Name of the vector variable
   * @param vector Vector field data
   * @return VTK unstructured grid object
   */
  vtkSmartPointer<vtkUnstructuredGrid> vector_to_vtk_unstructured_grid_data(const std::string& varname, const ndarray<F>&) const;
  // vtkSmartPointer<vtkUnstructuredGrid> scalars_to_vtk_unstructured_grid_data(
  //     const std::vector<std::string>& varname, const std::vector<ndarray<F>>& scalar) const;

  /**
   * @brief Convert mesh to VTK unstructured grid format
   * @return VTK unstructured grid representation of the mesh
   */
  virtual vtkSmartPointer<vtkUnstructuredGrid> to_vtu() const = 0;

  /**
   * @brief Initialize mesh from VTK unstructured grid
   * @param grid VTK unstructured grid to import
   */
  virtual void from_vtu(vtkSmartPointer<vtkUnstructuredGrid> grid) = 0;

  /**
   * @brief Check if a VTK grid is a simplicial mesh and get its dimension
   * @param grid VTK unstructured grid to check
   * @return 0 if non-simplicial, 2 for 2D triangular mesh, 3 for 3D tetrahedral mesh
   */
  static int check_simplicial_mesh_dims(vtkSmartPointer<vtkUnstructuredGrid> grid); // 0: nonsimplicial, 2 or 3: 2D or 3D
#endif
};

//// 
template <typename I, typename F>
void simplicial_unstructured_mesh<I, F>::from_legacy_vtk_file(const std::string& filename)
{
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkUnstructuredGridReader::New();
  reader->SetFileName(filename.c_str());
  reader->Update();
  from_vtu(reader->GetOutput());
#else
  ftk::fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

template <typename I, typename F>
void simplicial_unstructured_mesh<I, F>::from_vtk_unstructured_grid_file(const std::string& filename)
{
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkXMLUnstructuredGridReader::New();
  reader->SetFileName(filename.c_str());
  reader->Update();
  from_vtu(reader->GetOutput());
#else
  ftk::fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

template <typename I, typename F>
void simplicial_unstructured_mesh<I, F>::to_vtk_unstructured_grid_file(const std::string& filename) const
{
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( to_vtu() );
  writer->Write();
#else
  ftk::fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

template <typename I, typename F>
void simplicial_unstructured_mesh<I, F>::scalar_to_vtk_unstructured_grid_data_file(
    const std::string& filename, 
    const std::string& varname, 
    const ndarray<F>& scalar) const
{
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( scalar_to_vtk_unstructured_grid_data(varname, scalar) );
  writer->Write();
#else
  ftk::fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

template <typename I, typename F>
void simplicial_unstructured_mesh<I, F>::vector_to_vtk_unstructured_grid_data_file(
    const std::string& filename, 
    const std::string& varname, 
    const ndarray<F>& vector) const
{
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( vector_to_vtk_unstructured_grid_data(varname, vector) );
  writer->Write();
#else
  ftk::fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

#if FTK_HAVE_VTK
template <typename I, typename F>
vtkSmartPointer<vtkUnstructuredGrid> simplicial_unstructured_mesh<I, F>::scalar_to_vtk_unstructured_grid_data(
    const std::string& varname, const ndarray<F>& scalar) const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = to_vtu();

  vtkSmartPointer<vtkDataArray> array = vtkDoubleArray::New();
  array->SetName(varname.c_str());
  array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(n(0));
  for (int i = 0; i<n(0); i ++)
    array->SetTuple1(i, scalar[i]);

  grid->GetPointData()->AddArray(array);
  grid->GetPointData()->SetActiveScalars(varname.c_str());

  return grid;
}

template <typename I, typename F>
vtkSmartPointer<vtkUnstructuredGrid> simplicial_unstructured_mesh<I, F>::vector_to_vtk_unstructured_grid_data(
    const std::string& varname, const ndarray<F>& vector) const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = to_vtu();

  vtkSmartPointer<vtkDataArray> array = vtkDoubleArray::New();
  array->SetName(varname.c_str());
  array->SetNumberOfComponents(3); // TODO
  array->SetNumberOfTuples(n(0));

  if (vector.dim(0) == 2) {
    for (int i = 0; i<n(0); i ++)
      array->SetTuple3(i, vector(0, i), vector(1, i), 0);
  } else {
    for (int i = 0; i<n(0); i ++)
      array->SetTuple3(i, vector(0, i), vector(1, i), vector(2, i));
  }

  grid->GetPointData()->AddArray(array);
  grid->GetPointData()->SetActiveScalars(varname.c_str());

  return grid;
}

template <typename I, typename F>
int simplicial_unstructured_mesh<I, F>::check_simplicial_mesh_dims(vtkSmartPointer<vtkUnstructuredGrid> grid)
{
  int nd = 0;
#if (VTK_VERSION_NUMBER >= VTK_VERSION_CHECK(9, 2, 0))
  vtkSmartPointer<vtkUnsignedCharArray> types = grid->GetDistinctCellTypesArray();
  const int ntypes = types->GetNumberOfValues();

  if (ntypes == 1) {
    unsigned char type = types->GetValue(0);
    if (type == VTK_TRIANGLE) nd = 2;
    else if (type == VTK_TETRA) nd = 3;
    else nd = 0; // nonsimplicial
  } else 
    nd = 0; // nonsimplicial
#else 
  vtkSmartPointer<vtkCellTypes> types = vtkSmartPointer<vtkCellTypes>::New();
  grid->GetCellTypes(types);
  const int ntypes = types->GetNumberOfTypes();

  if (ntypes == 1) {
    unsigned char type = types->GetCellType(0);
    if (type == VTK_TRIANGLE) nd = 2;
    else if (type == VTK_TETRA) nd = 3;
    else nd = 0; // nonsimplicial
  } else 
    nd = 0; // nonsimplicial
#endif
  return nd;
}
#endif

} // namespace ftk

#endif
