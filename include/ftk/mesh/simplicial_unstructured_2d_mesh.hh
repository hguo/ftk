#ifndef _HYPERMESH_SIMPLEX_2D_HH
#define _HYPERMESH_SIMPLEX_2D_HH

#include <ftk/ftk_config.hh>
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
#include <vtkPointData.h>
#include <vtkPoints2D.h>
#endif

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_unstructured_2d_mesh { // 2D triangular mesh
  simplicial_unstructured_2d_mesh() {}

  simplicial_unstructured_2d_mesh(
      const std::vector<F>& coords, // coordinates of vertices; the dimension of the array is 2 * n_vertices
      const std::vector<I>& triangles); // vertex id of triangles; the dimension of the array is 3 * n_triangles

  simplicial_unstructured_2d_mesh(
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
  void smooth_scalar_gradient_jacobian(
      const ndarray<F>& f, 
      const F sigma,
      ndarray<F>& fs, // smoothed scalar field
      ndarray<F>& g,  // smoothed gradient field
      ndarray<F>& j   // smoothed jacobian field
  ) const; 

public: // io
  void scalar_to_vtk_unstructured_grid_data_file(const std::string& filename, const std::string& varname, const ndarray<F>&) const;
  void vector_to_vtk_unstructured_grid_data_file(const std::string& filename, const std::string& varname, const ndarray<F>&) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_vtk_unstructured_grid() const;
  vtkSmartPointer<vtkUnstructuredGrid> scalar_to_vtk_unstructured_grid_data(const std::string& varname, const ndarray<F>&) const;
  vtkSmartPointer<vtkUnstructuredGrid> vector_to_vtk_unstructured_grid_data(const std::string& varname, const ndarray<F>&) const;
  // vtkSmartPointer<vtkUnstructuredGrid> scalars_to_vtk_unstructured_grid_data(
  //     const std::vector<std::string>& varname, const std::vector<ndarray<F>>& scalar) const;
#endif
  void write_smoothing_kernel(const std::string& filename);
  bool read_smoothing_kernel(const std::string& filename);

public: // element iteration
  void element_for(int d, std::function<void(I)> f);

public: // mesh access
  std::set<I> sides(int d, I i);
  std::set<I> side_of(int d, I i);

  void get_triangle(I i, I tri[]) const;
  void get_edge(I i, I edge[]) const;
  void get_coords(I i, F coords[]) const;

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
simplicial_unstructured_2d_mesh<I, F>::simplicial_unstructured_2d_mesh(const std::vector<F> &coords_, const std::vector<I> &triangles_)
{
  vertex_coords.copy_vector(coords_);
  vertex_coords.reshape({2, coords_.size()/2});
  
  triangles.copy_vector(triangles_);
  triangles.reshape({3, triangles_.size()/3});
}

template <typename I, typename F>
size_t simplicial_unstructured_2d_mesh<I, F>::n(int d) const
{
  if (d == 0) return vertex_coords.dim(1);
  else if (d == 1) {
    return edges.dim(1);
  } else if (d == 2)
    return triangles.dim(1);
  else return 0;
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::build_vertex_links()
{
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::build_edges()
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
inline void simplicial_unstructured_2d_mesh<I, F>::build_smoothing_kernel(const F sigma)
{
  const F sigma2 = sigma * sigma;
  const F limit = F(3) * sigma;

  auto neighbors = [&](I i) { return vertex_edge_vertex[i]; };
  smoothing_kernel.resize(n(0));

  for (auto i = 0; i < n(0); i ++) {
    std::set<I> set;
    const F xi[2] = {vertex_coords[i*2], vertex_coords[i*2+1]};
    // fprintf(stderr, "i=%d, x={%f, %f}\n", i, xi[0], xi[1]);
    auto criteron = [&](I j) {
      const F xj[2] = {vertex_coords[j*2], vertex_coords[j*2+1]};
      if (vector_dist_2norm_2(xi, xj) < limit)
        return true;
      else return false;
    };

    auto operation = [&set](I j) {set.insert(j);};
    bfs<I, std::set<I>>(i, neighbors, operation, criteron);

    auto &kernel = smoothing_kernel[i];
    for (auto k : set) {
      const F xk[2] = {vertex_coords[2*k], vertex_coords[2*k+1]};
      const F d = vector_dist_2norm_2(xi, xk);
      const F w = std::exp(-(d*d) / (2*sigma*sigma)) / (sigma * std::sqrt(2.0 * M_PI));
      // fprintf(stderr, "d2=%f, w=%f\n", d2, w);
      kernel.push_back( std::make_tuple(k, w) );
      fprintf(stderr, "i=%d, k=%d, %f\n", i, k, w); 
    }

    // normalization
#if 0
    F sum = 0;
    for (int k = 0; k < kernel.size(); k ++)
      sum += std::get<1>(kernel[k]);
    for (int k = 0; k < kernel.size(); k ++) {
      std::get<1>(kernel[k]) /= sum;
      fprintf(stderr, "i=%d, k=%d, %f\n", i, k, std::get<1>(kernel[k]));// kernel.size());
    }
#endif
  }
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::smooth_scalar_gradient_jacobian(
    const ndarray<F>& f, const F sigma, 
    ndarray<F>& scalar, // smoothed scalar field
    ndarray<F>& grad,  // smoothed gradient field
    ndarray<F>& J) const // smoothed jacobian field
{
  const F sigma2 = sigma * sigma, 
          sigma4 = sigma2 * sigma2;

  scalar.reshape({n(0)});
  grad.reshape({2, n(0)});
  J.reshape({2, 2, n(0)});

  for (auto i = 0; i < smoothing_kernel.size(); i ++) {
    for (auto j = 0; j < smoothing_kernel[i].size(); j ++) {
      auto tuple = smoothing_kernel[i][j];
      const auto k = std::get<0>(tuple);
      const auto w = std::get<1>(tuple);
    
      const F d[2] = {vertex_coords[k*2] - vertex_coords[i*2], 
                      vertex_coords[k*2+1] - vertex_coords[i*2+1]};
      // const F r2 = d[0]*d[0] + d[1]*d[1];
      // const F r = std::sqrt(r2);

      // scalar
      scalar[i] += f[k] * w;

      // gradient
      grad(0, i) += - f[k] * w * d[0] / sigma2;
      grad(1, i) += - f[k] * w * d[1] / sigma2;

      // jacobian
      J(0, 0, i) += (d[0]*d[0] / sigma2 - 1) / sigma2 * f[k] * w;
      J(0, 1, i) += d[0]*d[1] / sigma4 * f[k] * w;
      J(1, 0, i) += d[0]*d[1] / sigma4 * f[k] * w;
      J(1, 1, i) += (d[1]*d[1] / sigma2 - 1) / sigma2 * f[k] * w;
    }
  }
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_2d_mesh<I, F>::sides(int d, I i)
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
std::set<I> simplicial_unstructured_2d_mesh<I, F>::side_of(int d, I i)
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
vtkSmartPointer<vtkUnstructuredGrid> simplicial_unstructured_2d_mesh<I, F>::to_vtk_unstructured_grid() const
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
  return grid;
}

template <typename I, typename F>
vtkSmartPointer<vtkUnstructuredGrid> simplicial_unstructured_2d_mesh<I, F>::scalar_to_vtk_unstructured_grid_data(
    const std::string& varname, const ndarray<F>& scalar) const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = to_vtk_unstructured_grid();

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
vtkSmartPointer<vtkUnstructuredGrid> simplicial_unstructured_2d_mesh<I, F>::vector_to_vtk_unstructured_grid_data(
    const std::string& varname, const ndarray<F>& vector) const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = to_vtk_unstructured_grid();

  vtkSmartPointer<vtkDataArray> array = vtkDoubleArray::New();
  array->SetName(varname.c_str());
  array->SetNumberOfComponents(3); // TODO
  array->SetNumberOfTuples(n(0));
  for (int i = 0; i<n(0); i ++) {
    array->SetTuple3(i, vector(0, i), vector(1, i), 0);
  }

  grid->GetPointData()->AddArray(array);
  grid->GetPointData()->SetActiveScalars(varname.c_str());

  return grid;
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::scalar_to_vtk_unstructured_grid_data_file(
    const std::string& filename, 
    const std::string& varname, 
    const ndarray<F>& scalar) const
{
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( scalar_to_vtk_unstructured_grid_data(varname, scalar) );
  writer->Write();
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::vector_to_vtk_unstructured_grid_data_file(
    const std::string& filename, 
    const std::string& varname, 
    const ndarray<F>& vector) const
{
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( vector_to_vtk_unstructured_grid_data(varname, vector) );
  writer->Write();
}
#else
template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::scalar_to_vtk_unstructured_grid_data_file(
    const std::string&, const std::string&, const ndarray<F>&) const
{
  fatal("FTK not compiled with VTK";
}
#endif

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::write_smoothing_kernel(const std::string& f)
{
  diy::serializeToFile(smoothing_kernel, f);
}

template <typename I, typename F>
bool simplicial_unstructured_2d_mesh<I, F>::read_smoothing_kernel(const std::string& f)
{
  return diy::unserializeFromFile(f, smoothing_kernel);
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::element_for(int d, std::function<void(I)> f)
{
  for (auto i = 0; i < n(d); i ++)
    f(i);
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::get_triangle(I i, I tri[]) const
{
  tri[0] = triangles(0, i);
  tri[1] = triangles(1, i);
  tri[2] = triangles(2, i);
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::get_edge(I i, I v[]) const
{
  v[0] = edges(0, i);
  v[1] = edges(1, i);
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::get_coords(I i, F coords[]) const
{
  coords[0] = vertex_coords(0, i);
  coords[1] = vertex_coords(1, i);
}

}
#endif
