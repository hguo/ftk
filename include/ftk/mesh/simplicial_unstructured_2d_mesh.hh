#ifndef _HYPERMESH_SIMPLEX_2D_HH
#define _HYPERMESH_SIMPLEX_2D_HH

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
#include <vtkUnstructuredGridReader.h>
#include <vtkPointData.h>
#include <vtkPoints2D.h>
#endif

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_unstructured_2d_mesh : public object { // 2D triangular mesh
  simplicial_unstructured_2d_mesh() {}

  simplicial_unstructured_2d_mesh(
      const std::vector<F>& coords, // coordinates of vertices; the dimension of the array is 2 * n_vertices
      const std::vector<I>& triangles); // vertex id of triangles; the dimension of the array is 3 * n_triangles

  simplicial_unstructured_2d_mesh(
      const ndarray<F>& coords_, // 2 * n_vertices
      const ndarray<I>& triangles_) // 3 * n_triangles
    : vertex_coords(coords_), triangles(triangles_) {build_triangles();}

  // dimensionality of the mesh
  int nd() const {return 2;}

  // numer of d-dimensional elements
  size_t n(int d) const;

  void build_edges();
  void build_triangles();

  void build_smoothing_kernel(F sigma);
  void smooth_scalar_gradient_jacobian(
      const ndarray<F>& f, 
      const F sigma,
      ndarray<F>& fs, // smoothed scalar field
      ndarray<F>& g,  // smoothed gradient field
      ndarray<F>& j   // smoothed jacobian field
  ) const; 

public: // io
  void from_vtk_file(const std::string& filename);

  void scalar_to_vtk_unstructured_grid_data_file(const std::string& filename, const std::string& varname, const ndarray<F>&) const;
  void vector_to_vtk_unstructured_grid_data_file(const std::string& filename, const std::string& varname, const ndarray<F>&) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_vtk_unstructured_grid() const;
  vtkSmartPointer<vtkUnstructuredGrid> scalar_to_vtk_unstructured_grid_data(const std::string& varname, const ndarray<F>&) const;
  vtkSmartPointer<vtkUnstructuredGrid> vector_to_vtk_unstructured_grid_data(const std::string& varname, const ndarray<F>&) const;
  // vtkSmartPointer<vtkUnstructuredGrid> scalars_to_vtk_unstructured_grid_data(
  //     const std::vector<std::string>& varname, const std::vector<ndarray<F>>& scalar) const;

  void from_vtk_unstructured_grid(vtkSmartPointer<vtkUnstructuredGrid> grid);
#endif
  void write_smoothing_kernel(const std::string& filename);
  bool read_smoothing_kernel(const std::string& filename);

public: // element iteration
  void element_for(int d, std::function<void(I)> f);

public: // mesh access
  std::set<I> sides(int d, I i) const;
  std::set<I> side_of(int d, I i) const;

  // std::set<I> side_of2(const I v[2]) const;

  void get_triangle(I i, I tri[]) const;
  void get_edge(I i, I edge[]) const;
  void get_coords(I i, F coords[]) const;

  bool find_edge(const I v[2], I& i) const;
  bool find_triangle(const I v[3], I& i) const;

private: // mesh connectivities
  ndarray<F> vertex_coords; // 2 * n_vertices
  std::vector<std::set<I>> vertex_side_of;
  std::vector<std::set<I>> vertex_edge_vertex;

  ndarray<I> edges; // 2 * n_edges
  ndarray<I> edges_side_of; // 2 * n_edges
  // std::map<std::tuple<I, I>, std::set<I>> edges_side_of;

  ndarray<I> triangles; // 3 * n_triangles
  ndarray<I> triangle_sides; // 3 * n_triangles

public: // additional mesh info
  std::map<std::tuple<I, I>, int> edge_id_map;
  std::map<std::tuple<I, I, I>, int> triangle_id_map;

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

  build_triangles();
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
bool simplicial_unstructured_2d_mesh<I, F>::find_edge(const I v[2], I &i) const
{
  const auto it = edge_id_map.find( std::make_tuple(v[0], v[1]) );
  if (it == edge_id_map.end()) return false;
  else {
    i = it->second;
    // fprintf(stderr, "edges(0,i)=%d, v[0]=%d\n", edges(0, i), v[0]);
    // fprintf(stderr, "edges(1,i)=%d, v[1]=%d\n", edges(1, i), v[1]);
    assert(edges(0, i) == v[0]);
    assert(edges(1, i) == v[1]);
    return true;
  }
}

template <typename I, typename F>
bool simplicial_unstructured_2d_mesh<I, F>::find_triangle(const I v[3], I &i) const
{
  const auto it = triangle_id_map.find( std::make_tuple(v[0], v[1], v[2]) );
  if (it == triangle_id_map.end()) return false;
  else {
    i = it->second;
    assert(triangles(0, i) == v[0]);
    assert(triangles(1, i) == v[1]);
    assert(triangles(2, i) == v[2]);
    return true;
  }
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::build_triangles()
{
  for (auto i = 0; i < n(2); i ++) {
    I v[3];
    get_triangle(i, v);
    std::sort(v, v+3);
    for (auto j = 0; j < 3; j ++)
      triangles(j, i) = v[j];

    triangle_id_map[ std::make_tuple(v[0], v[1], v[2]) ] = i;
  }
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::build_edges()
{
  // std::map<std::tuple<I, I>, std::set<I>> map_edges_side_of;
  std::map<I, std::set<I>> map_edges_side_of;

  int edge_count = 0;
  auto add_edge = [&](I tid, I v0, I v1) {
    if (v0 > v1) std::swap(v0, v1);
    const auto edge = std::make_tuple(v0, v1);
    I id;
    if (edge_id_map.find(edge) == edge_id_map.end()) {
      id = edge_count ++;
      edge_id_map[edge] = id;
    } else 
      id = edge_id_map[edge];

    map_edges_side_of[id].insert(tid);
      
    // fprintf(stderr, "adding edge #%d (%d, %d) from triangle %d\n", id, v0, v1, tid);
    // if (v0 == 0 && v1 == 1) {
    //   fprintf(stderr, "{0, 1}, tid=%d\n", tid);
    //   exit(1);
    // }
  };

  // fprintf(stderr, "triangles_88078=%d, %d, %d\n", triangles(0, 88078), triangles(1, 88078), triangles(2, 88078)); 
  for (auto i = 0; i < n(2); i ++) {
    add_edge(i, triangles(0, i), triangles(1, i));
    add_edge(i, triangles(1, i), triangles(2, i));
    add_edge(i, triangles(2, i), triangles(0, i));
  }

  edges.reshape(2, edge_id_map.size());
  edges_side_of.reshape({2, edge_id_map.size()}, -1);
  // int i = 0;
  for (const auto &kv : edge_id_map) {
    edges(0, kv.second) = std::get<0>(kv.first);
    edges(1, kv.second) = std::get<1>(kv.first);

    int j = 0;
    for (const auto tid : map_edges_side_of[kv.second])
    // for (const auto tid : map_edges_side_of[i])
      edges_side_of(j++, kv.second) = tid;
  }

  // for (auto i = 0; i < edges_side_of.dim(1); i ++)
  //   fprintf(stderr, "i=%d, tri0=%d, tri1=%d\n", i, edges_side_of(0, i), edges_side_of(1, i));

#if 0
  {
    I v[2] = {0, 1};
    I e;
    find_edge(v, e);
    auto s = side_of(1, e);
    // auto s = edges_side_of[std::make_tuple(v[0], v[1])];
    fprintf(stderr, "%zu, %d\n", s.size(), *s.begin());
  }
#endif
  // exit(1);

  vertex_side_of.resize(vertex_coords.dim(1));
  // fprintf(stderr, "resizing vertex_side_of, %zu\n", vertex_coords.dim(1));

  vertex_edge_vertex.resize(n(0));
  for (const auto &kv : edge_id_map) {
    const auto v0 = std::get<0>(kv.first),
               v1 = std::get<1>(kv.first);
    vertex_side_of[v0].insert(kv.second);
    vertex_side_of[v1].insert(kv.second);
  
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

  // for (auto i = 0; i < smoothing_kernel.size(); i ++) {
  parallel_for(smoothing_kernel.size(), std::thread::hardware_concurrency(), [&](int i) {
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
  });
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_2d_mesh<I, F>::sides(int d, I i) const
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
std::set<I> simplicial_unstructured_2d_mesh<I, F>::side_of(int d, I i) const
{
  if (d == 0)
    return vertex_side_of[i];
  else if (d == 1) {
    std::set<I> results;
    for (int j = 0; j < 2; j ++) 
      if (edges_side_of(j, i) >= 0) 
        results.insert(edges_side_of(j, i));
    return results;
#if 0
    int v[2];
    find_edge(v, i);
   
    const auto it = edges_side_of.find(std::make_tuple(v[0], v[1]));
    if (it == edges_side_of.end()) return std::set<I>();
    else return it->second;
#endif
  } else 
    return std::set<I>();
}

#if FTK_HAVE_VTK
template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::from_vtk_file(const std::string& filename)
{
  vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkUnstructuredGridReader::New();
  reader->SetFileName(filename.c_str());
  reader->Update();
  from_vtk_unstructured_grid(reader->GetOutput());
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::from_vtk_unstructured_grid(vtkSmartPointer<vtkUnstructuredGrid> grid)
{
  vtkIdType ncells = grid->GetNumberOfCells();
  std::vector<int> m_triangles;
  for (vtkIdType i = 0; i < ncells; i ++) {
    vtkSmartPointer<vtkCell> cell = grid->GetCell(i);
    if (cell->GetCellType() == VTK_TRIANGLE) {
      vtkIdType v[3] = {cell->GetPointId(0), cell->GetPointId(1), cell->GetPointId(2)};
      std::sort(v, v+3);
      for (int j = 0; j < 3; j ++)
        m_triangles.push_back(v[j]);
    }
  }
  triangles.reshape({3, m_triangles.size()/3});
  triangles.from_vector(m_triangles);

  vtkIdType npts = grid->GetNumberOfPoints();
  vertex_coords.reshape({2, size_t(npts)});
  for (vtkIdType i = 0; i < npts; i ++) {
    double x[3];
    grid->GetPoint(i, x);
    vertex_coords(0, i) = x[0];
    vertex_coords(1, i) = x[1];
  }

  build_triangles();
  build_edges();
}

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
  fatal("FTK not compiled with VTK");
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::from_vtk_file(const std::string& filename)
{
  fatal("FTK not compiled with VTK");
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

#if 0
template <typename I, typename F>
std::set<I> simplicial_unstructured_2d_mesh<I, F>::side_of2(const I v[2]) const
{
  const auto edge = std::make_tuple(v[0], v[1]);
  const auto it = edges_side_of.find(edge);
  if (it == edges_side_of.end()) return std::set<I>();
  else return it->second;
}
#endif

}
#endif
