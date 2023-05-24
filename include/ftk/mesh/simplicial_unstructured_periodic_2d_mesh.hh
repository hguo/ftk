#ifndef _FTK_MESH_UNSTRUCTURED_PERIODIC_2D_HH
#define _FTK_MESH_UNSTRUCTURED_PERIODIC_2D_HH

#include <ftk/config.hh>
#include <ftk/mesh/simplicial_unstructured_mesh.hh>
#include <ftk/mesh/simplicial_unstructured_3d_mesh.hh>
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>
#include <ftk/utils/string.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_unstructured_periodic_2d_mesh : public object {
  simplicial_unstructured_periodic_2d_mesh(
      std::shared_ptr<simplicial_unstructured_2d_mesh<I, F>> m2_, 
      std::shared_ptr<simplicial_unstructured_3d_mesh<I, F>> m3_) : m2(m2_), m3(m3_) { initialize(); }

  void initialize();

  size_t n(int d) const;
  size_t n_ordinal(int d) const;
  size_t n_interval(int d) const;
  
  I flat_vertex_id(I i) const { return mod(i, m2->n(0)); }
  I flat_vertex_time(I i) const { return i / m2->n(0); }
  
  void get_simplex(int d, I i, I v[]) const;
  void get_coords(I i, F coords[]) const;
  
  std::set<I> sides(int d, I i) const;
  std::set<I> side_of(int d, I i) const;
  
  void element_for(int d, std::function<void(I)> f, 
      int xl = FTK_XL_NONE, int nthreads=std::thread::hardware_concurrency(), bool affinity = false) const;
  
  void element_for_ordinal(int d, int t, std::function<void(I)> f,
      int xl = FTK_XL_NONE, int nthreads=std::thread::hardware_concurrency(), bool affinity = false) const;
  
  void element_for_interval(int d, int t, std::function<void(I)> f,
      int xl = FTK_XL_NONE, int nthreads=std::thread::hardware_concurrency(), bool affinity = false) const;

public:
  void to_vtu_file(const std::string& filename, int np) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_vtu(int np) const;
#endif
  
  template <typename T>
  void array_to_vtu(const std::string& filename, 
      const std::vector<std::string>& varnames, 
      const std::vector<ndarray<T>>& arrays) const;

protected:
  I normalize(int d, I k) const;
  
  static I mod(I val, I m);

protected:
  // the base mesh with identical lower and upper 2D meshes
  std::shared_ptr<simplicial_unstructured_3d_mesh<I, F>> m3; 
  std::shared_ptr<simplicial_unstructured_2d_mesh<I, F>> m2;

  // std::vector<I> m3_triangles, m3_edges;
  std::vector<I> m3_ordinal_triangles, m3_ordinal_edges;
  std::vector<I> m3_interval_triangles, m3_interval_edges;
};

/////////
template <typename I, typename F> 
size_t simplicial_unstructured_periodic_2d_mesh<I, F>::n(int d) const // number of elements per unit period
{
  if (d == 3) 
    return m3->n(3);
  if (d == 2)
    return m3_ordinal_triangles.size() + m3_interval_triangles.size();
  else if (d == 1) 
    return m3_ordinal_edges.size() + m3_interval_edges.size();
  else // TODO
    return m3->n(d) - m2->n(d);
}

template <typename I, typename F> 
size_t simplicial_unstructured_periodic_2d_mesh<I, F>::n_ordinal(int d) const
{
  if (d == 3) return m3->n(d);
  else if (d == 2) return m3_ordinal_triangles.size();
  else if (d == 1) return m3_ordinal_edges.size();
  else return m2->n(d);
}

template <typename I, typename F> 
size_t simplicial_unstructured_periodic_2d_mesh<I, F>::n_interval(int d) const
{
  if (d == 2) return m3_interval_triangles.size();
  else return m3->n(d) - 2 * m2->n(d);
}

template <typename I, typename F>
I simplicial_unstructured_periodic_2d_mesh<I, F>::mod(I v, I m)
{
  I mod = v % m;
  if (v < 0) mod += m;
  return mod;
}

template <typename I, typename F>
void simplicial_unstructured_periodic_2d_mesh<I, F>::element_for(int d, std::function<void(I)> f, 
    int xl, int nthreads, bool affinity) const
{
  parallel_for(n(d), [&](int i) {f(i);}, xl, nthreads, affinity);
  // for (auto i = 0; i < n(d); i ++)
  //   f(i);
}

template <typename I, typename F>
void simplicial_unstructured_periodic_2d_mesh<I, F>::element_for_ordinal(
    int d, int t, std::function<void(I)> f,
    int xl, int nthreads, bool affinity) const
{
  parallel_for(n_ordinal(d), [&](int i) {
        const auto id = i + t * n(d);
        f(id);
      }, xl, nthreads, affinity);
}

template <typename I, typename F>
void simplicial_unstructured_periodic_2d_mesh<I, F>::element_for_interval(int d, int t, std::function<void(I)> f,
    int xl, int nthreads, bool affinity) const
{
  parallel_for(n_interval(d), [&](int i) {
    const auto id = i + n_ordinal(d) + t * n(d);
    f(id);
  }, xl, nthreads, affinity);
}

template <typename I, typename F>
void simplicial_unstructured_periodic_2d_mesh<I, F>::get_simplex(int d, I k, I verts[]) const
{
  const I i = mod(k, n(d)), t = std::floor(double(k) / n(d));
  const I offset = t * n(0); // 0d gid offset for verts

  if (d == 3) {
    m3->get_simplex(d, i, verts);
  } else if (d == 2) {
    if (i < n_ordinal(d)) 
      m3->get_simplex(d, m3_ordinal_triangles[i], verts);
    else 
      m3->get_simplex(d, m3_interval_triangles[i - n_ordinal(d)], verts);
  } else if (d == 1) {
    if (i < n_ordinal(d)) 
      m3->get_simplex(d, m3_ordinal_edges[i], verts);
    else 
      m3->get_simplex(d, m3_interval_edges[i - n_ordinal(d)], verts);
  } else 
    assert(false); // not implemented yet
    
  for (int i = 0; i < d+1; i ++)
    verts[i] += offset;
}

template <typename I, typename F>
void simplicial_unstructured_periodic_2d_mesh<I, F>::get_coords(I i, F coords[]) const
{
  const I k = flat_vertex_id(i),
          t = flat_vertex_time(i);
  m2->get_coords(k, coords);

  coords[2] = t;
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_periodic_2d_mesh<I, F>::sides(int d, I k) const 
{
  const I i = mod(k, n(d)), t = std::floor((double)k / n(d));
  const I offset = t * n(d);
  std::set<I> results;

  std::set<I> m3sides = m3->sides(d, i);
  for (const auto side : m3sides)
    results.insert(normalize(d, side + offset));

  return results;
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_periodic_2d_mesh<I, F>::side_of(int d, I k) const
{
  const I i = mod(k, n(d)), t = std::floor((double)k / n(d));
  std::set<I> results;

  std::set<I> m3sideof = m3->side_of(d, i);
  for (const auto sideof : m3sideof)
    results.insert(normalize(d+1, sideof + t * this->n(d)));

  // check if the queried simplex is purely spatial
  // - if pure, we also need to find sideofs in the previous period
  // - otherwise, just return simplices in the current period (plus offset)
  bool pure = true;
  I verts[4];
  m3->get_simplex(d, i, verts);
  for (int i = 0; i < d; i ++)
    if (verts[i] >= m3->n(0)) // any of the verts is in the upper plane
      return pure = false;

  if (!pure) {
    std::set<I> m3sideof = m3->side_of(d, i + n(d));
    for (const auto sideof : m3sideof)
      results.insert(normalize(d+1, sideof + (t - 1) * this->n(d)));
  }
}

template <typename I, typename F>
I simplicial_unstructured_periodic_2d_mesh<I, F>::normalize(int d, I k) const
{
  if (d == 0 || d == 3) return k; // no need to normalize
  else { // edges or triangles
    I verts[4];
    const I i = mod(k, n(d)), t = std::floor(double(k) / n(d));
 
    bool upper = true;
    m3->get_simplex(d, i, verts);
    for (int i = 0; i < d; i ++) {
      if (verts[i] < m3->n(0)) { // any of the verts is in the lower plane
        upper = false;
      }
    }

    if (upper) {
      for (int i = 0; i < d; i ++)
        verts[i] = verts[i] - m3->n(0);

      I j;
      bool succ = m3->find_simplex(d, verts, j);
      assert(succ);
      
      return j + t * n(d);
    } else 
      return k;
  }

}

template <typename I, typename F>
void simplicial_unstructured_periodic_2d_mesh<I, F>::initialize()
{
  const I n0 = m2->n(0);
  const I m2n2 = m2->n(2);

  // extract all interval triangles and edges in m3
  for (auto i = 0; i < m3->n(2); i ++) {
    I verts[3];
    m3->get_simplex(2, i, verts);

    bool all_lower = true, all_upper = true;
    for (auto j = 0; j < 3; j ++) {
      if (verts[j] < n0) all_upper = false;
      else /* if (verts[j] >= n0) */ all_lower = false;
    }

    if (all_lower) {
      // fprintf(stderr, "adding lower, %zu, %zu, %zu\n", verts[0], verts[1], verts[2]);
      m3_ordinal_triangles.push_back(i);
    }
    else if (!all_upper)
      m3_interval_triangles.push_back(i);
  }

  // interval edges
  for (auto i = 0; i < m3->n(1); i ++) {
    I verts[2];
    m3->get_simplex(1, i, verts);
    
    bool all_lower = true, all_upper = true;
    for (auto j = 0; j < 2; j ++) {
      if (verts[j] < n0) all_upper = false;
      else /* if (verts[j] >= n0) */ all_lower = false;
    }

    if (all_lower) 
      m3_ordinal_edges.push_back(i);
    else if (!all_upper)
      m3_interval_edges.push_back(i);
  }

  // fprintf(stderr, "%zu, %zu\n",
  //     m3_ordinal_triangles.size(), m3_interval_triangles.size());
}

template <typename I, typename F>
void simplicial_unstructured_periodic_2d_mesh<I, F>::to_vtu_file(const std::string& filename, int np) const
{
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( to_vtu(np) );
  writer->Write();
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

#if FTK_HAVE_VTK
template <typename I, typename F>
vtkSmartPointer<vtkUnstructuredGrid> simplicial_unstructured_periodic_2d_mesh<I, F>::
to_vtu(int np) const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> pts = vtkPoints::New();
  pts->SetNumberOfPoints(n(0) * np);

  fprintf(stderr, "converting to vtu, n0=%zu, n3=%zu\n", n(0), n(3));

  for (int p = 0; p < np; p ++) {
    const I offset = p * n(0);
    for (I i=0; i < n(0); i++) {
      // pts->SetPoint(i, vertex_coords[i*3], vertex_coords[i*3+1], vertex_coords[i*3+2]); 
      F coords[3];
      get_coords(i + offset, coords);
      pts->SetPoint(i + offset, coords[0], coords[1], coords[2]);
      // fprintf(stderr, "%f, %f, %f\n", coords[0], coords[1], coords[2]);
    }
  }
  grid->SetPoints(pts);

  for (int p = 0; p < np-1; p ++) {
    const I offset = p * n(0);
    for (int i=0; i < n(3); i ++) {
      // vtkIdType ids[4] = {tetrahedra[i*4], tetrahedra[i*4+1], tetrahedra[i*4+2], tetrahedra[i*4+3]};
      I tet[4];
      m3->get_simplex(3, i, tet);
      vtkIdType ids[4] = {
        tet[0] + offset, 
        tet[1] + offset, 
        tet[2] + offset, 
        tet[3] + offset};
      // fprintf(stderr, "%d, %d, %d, %d\n", ids[0], ids[1], ids[2], ids[3]);
      // fprintf(stderr, "adding tet %d: %d, %d, %d, %d\n", i, tet[0], tet[1], tet[2], tet[3]);
      grid->InsertNextCell(VTK_TETRA, 4, ids);
    }
  }

  // grid->PrintSelf(std::cerr, vtkIndent(2));

  return grid;
}
#endif

template <typename I, typename F>
template <typename T>
void simplicial_unstructured_periodic_2d_mesh<I, F>::array_to_vtu(
    const std::string& filename, 
    const std::vector<std::string>& varnames, 
    const std::vector<ndarray<T>>& arrays) const
{
#if FTK_HAVE_VTK
  auto grid = to_vtu(arrays[0].dim(1));

  for (int i = 0; i < varnames.size(); i ++) {
    auto data = arrays[i].to_vtk_data_array(varnames[i]);
    grid->GetPointData()->AddArray(data);
  }

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInputData( grid );
  writer->Write();
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

}

#endif
