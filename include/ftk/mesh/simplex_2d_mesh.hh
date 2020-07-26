#ifndef _HYPERMESH_SIMPLEX_2D_HH
#define _HYPERMESH_SIMPLEX_2D_HH

#include <ftk/ftk_config.hh>
#include <ftk/ndarray.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/algorithms/bfs.hh>
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
#include <vtkPointData.h>
#include <vtkPoints2D.h>
#endif

namespace ftk {

template <typename I=int, typename F=double>
struct simplex_2d_mesh { // 2D triangular mesh
  simplex_2d_mesh() {}

  simplex_2d_mesh(
      const std::vector<F>& coords, // coordinates of vertices; the dimension of the array is 2 * n_vertices
      const std::vector<I>& triangles); // vertex id of triangles; the dimension of the array is 3 * n_triangles

  simplex_2d_mesh(
      const ndarray<F>& coords_, // 2 * n_vertices
      const ndarray<I>& triangles_) // 3 * n_triangles
    : vertex_coords(coords_), triangles(triangles_) {}

  // dimensionality of the mesh
  int nd() const {return 2;}

  // numer of d-dimensional elements
  size_t n(int d) const;

  void build_vertex_links();
  void build_edges();

  void build_smoothing_kernel(F sigma);
  ndarray<F> smooth_scalar_field(const ndarray<F> &f);

public: // io
  void scalar_to_vtk_unstructured_grid_data_file(const std::string& filename, const std::string& varname, const ndarray<F>& scalar) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> scalar_to_vtk_unstructured_grid_data(const std::string& varname, const ndarray<F>& scalar) const;
#endif

public: // mesh access
  std::set<I> sides(int d, I i);
  std::set<I> side_of(int d, I i);

private: // mesh connectivities
  ndarray<F> vertex_coords; // 2 * n_vertices
  std::vector<std::set<I>> vertex_side_of;
  std::vector<std::set<I>> vertex_edge_vertex;

  ndarray<I> edges; // 2 * n_edges
  ndarray<I> edges_side_of; // 2 * n_edges

  ndarray<I> triangles; // 3 * n_triangles
  ndarray<I> triangle_sides; // 3 * n_triangles

private:
  std::vector<std::vector<std::tuple<I/*vert*/, F/*weight*/>>> smoothing_kernel;
};

/////////

template <typename I, typename F>
simplex_2d_mesh<I, F>::simplex_2d_mesh(const std::vector<F> &coords_, const std::vector<I> &triangles_)
{
  vertex_coords.copy_vector(coords_);
  vertex_coords.reshape({2, coords_.size()/2});
  
  triangles.copy_vector(triangles_);
  triangles.reshape({3, triangles_.size()/3});
}

template <typename I, typename F>
size_t simplex_2d_mesh<I, F>::n(int d) const
{
  if (d == 0) return vertex_coords.dim(1);
  else if (d == 1) {
    return edges.dim(1);
  } else if (d == 2)
    return triangles.dim(1);
  else return 0;
}

template <typename I, typename F>
void simplex_2d_mesh<I, F>::build_vertex_links()
{
}

template <typename I, typename F>
void simplex_2d_mesh<I, F>::build_edges()
{
  typedef std::tuple<I, I> edge_t;
  std::set<edge_t> unique_edges;

  auto convert_edge = [](edge_t e) {
    if (std::get<0>(e) > std::get<1>(e))
      return std::make_tuple(std::get<1>(e), std::get<0>(e));
    else return e;
  };

  for (auto i = 0; i < triangles.dim(1); i ++) {
    unique_edges.insert(std::make_tuple(triangles(0, i), triangles(1, i)));
    unique_edges.insert(std::make_tuple(triangles(1, i), triangles(2, i)));
    unique_edges.insert(std::make_tuple(triangles(2, i), triangles(3, i)));
  }

  edges.reshape(2, unique_edges.size());
  vertex_side_of.resize(vertex_coords.dim(1));
  // fprintf(stderr, "resizing vertex_side_of, %zu\n", vertex_coords.dim(1));

  int i = 0;
  vertex_edge_vertex.resize(n(0));
  for (const auto e : unique_edges) {
    auto v0 = edges(0, i) = std::get<0>(e);
    auto v1 = edges(1, i) = std::get<1>(e);
    vertex_side_of[v0].insert(i);
    vertex_side_of[v1].insert(i);
    i ++;
  
    vertex_edge_vertex[v0].insert(v1);
    vertex_edge_vertex[v1].insert(v0);
  }
}

template <typename I, typename F>
inline void simplex_2d_mesh<I, F>::build_smoothing_kernel(const F sigma)
{
  const F limit = F(3) * sigma;

  auto neighbors = [&](I i) {
    return vertex_edge_vertex[i];
#if 0
    std::set<I> results;
    for (auto edge : side_of(0, i))
      for (auto side : sides(1, edge))
        results.insert(side);
    return results;
#endif
  };

  smoothing_kernel.resize(n(0));

  for (auto i = 0; i < n(0); i ++) {
    std::set<I> set;
    const F xi[2] = {vertex_coords(0, i), vertex_coords(1, i)};
    // fprintf(stderr, "i=%d, x={%f, %f}\n", i, xi[0], xi[1]);
    auto criteron = [&](I j) {
      const F xj[2] = {vertex_coords(0, j), vertex_coords(1, j)};
      if (vector_dist_2norm_2(xi, xj) < limit)
        return true;
      else return false;
    };

    auto operation = [&set](I j) {set.insert(j);};
    bfs<I, std::set<I>>(i, neighbors, operation, criteron);

    auto &kernel = smoothing_kernel[i];
    for (auto k : set) {
      const F xk[2] = {vertex_coords(0, k), vertex_coords(1, k)};
      const F d = vector_dist_2norm_2(xi, xk);
      const F w = std::exp(-(d*d) / (sigma*sigma)) / (sigma * std::sqrt(2.0 * M_PI));
      // fprintf(stderr, "d2=%f, w=%f\n", d2, w);
      kernel.push_back( std::make_tuple(k, w) );
    }

    // normalization
    F sum = 0;
    for (int k = 0; k < kernel.size(); k ++)
      sum += std::get<1>(kernel[k]);
    for (int k = 0; k < kernel.size(); k ++) {
      std::get<1>(kernel[k]) /= sum;
      fprintf(stderr, "i=%d, k=%d, %f\n", i, k, std::get<1>(kernel[k]));// kernel.size());
    }
  }
}

template <typename I, typename F>
ndarray<F> simplex_2d_mesh<I, F>::smooth_scalar_field(const ndarray<F>& f)
{
  ndarray<F> result;
  result.reshape(f);

  for (auto i = 0; i < smoothing_kernel.size(); i ++) {
    for (auto j = 0; j < smoothing_kernel[i].size(); j ++) {
      auto tuple = smoothing_kernel[i][j];
      const auto k = std::get<0>(tuple);
      const auto w = std::get<1>(tuple);
      result[i] += f[k] * w;
    }
  }

  return result;
}

template <typename I, typename F>
std::set<I> simplex_2d_mesh<I, F>::sides(int d, I i)
{
  std::set<I> results;
  if (d == 1) {
    results.insert( edges(0, i) );
    results.insert( edges(1, i) );
  } else if (d == 2) {
    results.insert( triangle_sides(0, i) );
    results.insert( triangle_sides(1, i) );
    results.insert( triangle_sides(2, i) );
  }
  return results;
}

template <typename I, typename F>
std::set<I> simplex_2d_mesh<I, F>::side_of(int d, I i)
{
  std::set<I> results;
  if (d == 0)
    return vertex_side_of[i];
  else if (d == 1) {
    results.insert(edges_side_of(0, i));
    results.insert(edges_side_of(1, i));
  }
  return results;
}

#if FTK_HAVE_VTK
template <typename I, typename F>
vtkSmartPointer<vtkUnstructuredGrid> simplex_2d_mesh<I, F>::scalar_to_vtk_unstructured_grid_data(
    const std::string& varname, const ndarray<F>& scalar) const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> pts = vtkPoints::New();
  pts->SetNumberOfPoints(n(0));

  for (int i=0; i<n(0); i++)
    pts->SetPoint(i, vertex_coords[i*2], vertex_coords[i*2+1], 0); 

  for (int i=0; i<n(2); i ++) {
    vtkIdType ids[3] = {triangles[i*3], triangles[i*3+1], triangles[i*3+2]};
    grid->InsertNextCell(VTK_TRIANGLE, 3, ids);
  }

  grid->SetPoints(pts);

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
void simplex_2d_mesh<I, F>::scalar_to_vtk_unstructured_grid_data_file(
    const std::string& filename, 
    const std::string& varname, 
    const ndarray<F>& scalar) const
{
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( scalar_to_vtk_unstructured_grid_data(varname, scalar) );
  writer->Write();
}
#endif

}

#endif
