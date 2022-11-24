#ifndef _FTK_MESH_UNSTRUCTURED_3D_HH
#define _FTK_MESH_UNSTRUCTURED_3D_HH

#include <ftk/config.hh>
#include <ftk/mesh/simplicial_unstructured_mesh.hh>
#include <ftk/mesh/simplicial_unstructured_extruded_2d_mesh.hh>
#include <ftk/utils/string.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_unstructured_3d_mesh : public simplicial_unstructured_mesh<I, F> {
  friend class diy::Serialization<simplicial_unstructured_3d_mesh<I, F>>;
  
  simplicial_unstructured_3d_mesh() {}

  simplicial_unstructured_3d_mesh(
      const std::vector<F>& coords, 
      const std::vector<I>& tetrahedra);

  int nd() const {return 3;}
  virtual size_t n(int d, bool part = false) const;
  
  void build_smoothing_kernel(F sigma);
  void smooth_scalar_gradient_jacobian(
      const ndarray<F>& f, 
      const F sigma,
      ndarray<F>& fs, // smoothed scalar field
      ndarray<F>& g,  // smoothed gradient field
      ndarray<F>& j   // smoothed jacobian field
  ) const; 

public: // io
  static std::shared_ptr<simplicial_unstructured_3d_mesh<I, F>> from_file(const std::string& filename);
  void to_file(const std::string &filename) const;

  void to_vtu_file(const std::string& filename) const;
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_vtu() const;
  void from_vtu(vtkSmartPointer<vtkUnstructuredGrid> grid);
#endif

public: 
  virtual void element_for(int d, std::function<void(I)> f) {} // TODO

private: // use get_simplex() and find_simplex instead
  void get_tetrahedron(I i, I tet[]) const;
  void get_triangle(I i, I tri[]) const;
  void get_edge(I i, I edge[]) const;

  bool find_tetrahedron(const I v[4], I& i) const;
  bool find_triangle(const I v[3], I& i) const;
  bool find_edge(const I v[2], I& i) const;

public:
  virtual std::set<I> sides(int d, I i) const { return std::set<int>(); } // TODO
  virtual std::set<I> side_of(int d, I i) const;

  virtual void get_simplex(int d, I i, I simplex[]) const;
  virtual bool find_simplex(int d, const I verts[], I& i) const;
  virtual void get_coords(I i, F coords[]) const;

  virtual const ndarray<F>& get_coords() const {return vertex_coords;}

private:
  void initialize();

  void build_edges();
  void build_triangles();
  void build_tetrahedra();

private: // mesh conn
  ndarray<F> vertex_coords; // 3*n_vert
  std::vector<std::set<I>> vertex_side_of;
  std::vector<std::set<I>> vertex_edge_vertex;

  ndarray<I> edges;
  std::vector<std::set<I>> edges_side_of;

  ndarray<I> triangles;
  ndarray<I> triangle_sides;
  // ndarray<I> triangle_side_of;
  std::vector<std::set<I>> triangle_side_of;

  ndarray<I> tetrahedra;
  ndarray<I> tetrahedra_sides;

public:
  std::map<std::tuple<I, I>, int> edge_id_map;
  std::map<std::tuple<I, I, I>, int> triangle_id_map;
  std::map<std::tuple<I, I, I, I>, int> tetrahedron_id_map;
  
  std::vector<std::vector<std::tuple<I/*vert*/, F/*weight*/>>> smoothing_kernel;
};

//////////
template <typename I, typename F>
simplicial_unstructured_3d_mesh<I, F>::simplicial_unstructured_3d_mesh(
    const std::vector<F>& coords_,
    const std::vector<I>& tetrahedra_)
{
  vertex_coords.copy_vector(coords_);
  vertex_coords.reshape({3, coords_.size()/3});
  
  tetrahedra.copy_vector(tetrahedra_);
  tetrahedra.reshape({4, tetrahedra_.size()/4});
}

template <typename I, typename F>
void simplicial_unstructured_3d_mesh<I, F>::build_tetrahedra()
{
  for (auto i = 0; i < n(3); i ++) {
    I v[4];
    get_tetrahedron(i, v);
    std::sort(v, v+4);
    for (auto j = 0; j < 4; j ++)
      tetrahedra(j, i) = v[j];

    tetrahedron_id_map[ std::make_tuple(v[0], v[1], v[2], v[3]) ] = i;
  }
}

template <typename I, typename F>
void simplicial_unstructured_3d_mesh<I, F>::build_triangles()
{
  std::map<I, std::set<I>> map_triangle_side_of;

  int triangle_count = 0;
  auto add_triangle = [&](I tid, I v0, I v1, I v2) {
    const auto tri = std::make_tuple(v0, v1, v2);
    I id;
    if (triangle_id_map.find(tri) == triangle_id_map.end()) {
      id = triangle_count ++;
      triangle_id_map[tri] = id;
    } else 
      id = triangle_id_map[tri];

    map_triangle_side_of[id].insert(tid);
  };

  for (auto i = 0; i < n(3); i ++) {
    add_triangle(i, tetrahedra(0, i), tetrahedra(1, i), tetrahedra(2, i));
    add_triangle(i, tetrahedra(0, i), tetrahedra(1, i), tetrahedra(3, i));
    add_triangle(i, tetrahedra(0, i), tetrahedra(2, i), tetrahedra(3, i));
    add_triangle(i, tetrahedra(1, i), tetrahedra(2, i), tetrahedra(3, i));
  }

  triangles.reshape(3, triangle_id_map.size());
  triangle_side_of.resize(triangle_id_map.size());
  for (const auto &kv : triangle_id_map) {
    triangles(0, kv.second) = std::get<0>(kv.first);
    triangles(1, kv.second) = std::get<1>(kv.first);
    triangles(2, kv.second) = std::get<2>(kv.first);

    for (const auto tid : map_triangle_side_of[kv.second])
      triangle_side_of[kv.second].insert(tid);
  }
}

template <typename I, typename F>
void simplicial_unstructured_3d_mesh<I, F>::build_edges()
{
  std::map<I, std::set<I>> map_edges_side_of;

  int edge_count = 0;
  auto add_edge = [&](I tid, I v0, I v1) {
    const auto edge = std::make_tuple(v0, v1);
    I id;
    if (edge_id_map.find(edge) == edge_id_map.end()) {
      id = edge_count ++;
      edge_id_map[edge] = id;
    } else 
      id = edge_id_map[edge];

    map_edges_side_of[id].insert(tid);
  };

  // fprintf(stderr, "triangles_88078=%d, %d, %d\n", triangles(0, 88078), triangles(1, 88078), triangles(2, 88078)); 
  for (auto i = 0; i < n(2); i ++) {
    add_edge(i, triangles(0, i), triangles(1, i));
    add_edge(i, triangles(0, i), triangles(2, i));
    add_edge(i, triangles(1, i), triangles(2, i));
  }

  edges.reshape(2, edge_id_map.size());
  edges_side_of.resize(edge_id_map.size());
  for (const auto &kv : edge_id_map) {
    edges(0, kv.second) = std::get<0>(kv.first);
    edges(1, kv.second) = std::get<1>(kv.first);

    for (const auto tid : map_edges_side_of[kv.second])
      edges_side_of[kv.second].insert(tid);
  }

  /////// 
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
size_t simplicial_unstructured_3d_mesh<I, F>::n(int d, bool part /* TODO */) const
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

template <typename I, typename F>
void simplicial_unstructured_3d_mesh<I, F>::initialize()
{
  build_tetrahedra();
  build_triangles();
  build_edges();
  
  fprintf(stderr, "3d mesh initialized: #tet=%zu, #tri=%zu, #edge=%zu, #vert=%zu\n", 
      tetrahedron_id_map.size(),
      triangle_id_map.size(), 
      edge_id_map.size(),
      vertex_edge_vertex.size());
}

template <typename I, typename F>
void simplicial_unstructured_3d_mesh<I, F>::to_vtu_file(const std::string& filename) const
{
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( to_vtu() );
  writer->Write();
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

#if FTK_HAVE_VTK  
template <typename I, typename F>
void simplicial_unstructured_3d_mesh<I, F>::from_vtu(vtkSmartPointer<vtkUnstructuredGrid> grid)
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

  initialize();
}

template <typename I, typename F>
vtkSmartPointer<vtkUnstructuredGrid> simplicial_unstructured_3d_mesh<I, F>::
to_vtu() const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> pts = vtkPoints::New();
  pts->SetNumberOfPoints(n(0));

  // fprintf(stderr, "converting to vtu, n0=%zu, n3=%zu\n", n(0), n(3));

  for (I i=0; i < n(0); i++) {
    // pts->SetPoint(i, vertex_coords[i*3], vertex_coords[i*3+1], vertex_coords[i*3+2]); 
    F coords[3];
    get_coords(i, coords);
    pts->SetPoint(i, coords[0], coords[1], coords[2]);
  }

  for (int i=0; i < n(3); i ++) {
    // vtkIdType ids[4] = {tetrahedra[i*4], tetrahedra[i*4+1], tetrahedra[i*4+2], tetrahedra[i*4+3]};
    I tet[4];
    get_simplex(3, i, tet);
    vtkIdType ids[4] = {tet[0], tet[1], tet[2], tet[3]};
    // fprintf(stderr, "adding tet %d: %d, %d, %d, %d\n", i, tet[0], tet[1], tet[2], tet[3]);
    grid->InsertNextCell(VTK_TETRA, 4, ids);
  }

  grid->SetPoints(pts);
  // grid->PrintSelf(std::cerr, vtkIndent(2));

  return grid;
}
#endif

template <typename I, typename F>
std::shared_ptr<simplicial_unstructured_3d_mesh<I, F>> simplicial_unstructured_3d_mesh<I, F>::from_file(const std::string& filename)
{
  std::shared_ptr<simplicial_unstructured_3d_mesh<I, F>> m(new simplicial_unstructured_3d_mesh<>);
  if (ends_with_lower(filename, "vtu")) {
#if FTK_HAVE_VTK
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkXMLUnstructuredGridReader::New();
    reader->SetFileName( filename.c_str() );
    reader->Update();
    m->from_vtu( reader->GetOutput() );
#else
    fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
  } else {
    diy::unserializeFromFile(filename, *m);
  }
  return m;
}

template <typename I, typename F>
void simplicial_unstructured_3d_mesh<I, F>::to_file(const std::string& filename) const
{
  if (ends_with_lower(filename, "vtu")) {
    to_vtu_file(filename);
  } else {
    diy::serializeToFile(*this, filename);
  }
}

template <typename I, typename F>
void simplicial_unstructured_3d_mesh<I, F>::get_simplex(int d, I i, I v[]) const
{
  if (d == 3) get_tetrahedron(i, v);
  else if (d == 2) get_triangle(i, v);
  else if (d == 1) get_edge(i, v);
}

template <typename I, typename F>
bool simplicial_unstructured_3d_mesh<I, F>::find_simplex(int d, const I v[], I &i) const
{
  if (d == 3) return find_tetrahedron(v, i);
  else if (d == 2) return find_triangle(v, i);
  else if (d == 1) return find_edge(v, i);
  else return false;
}

template <typename I, typename F>
void simplicial_unstructured_3d_mesh<I, F>::get_tetrahedron(I i, I tet[]) const
{
  for (int j = 0; j < 4; j ++)
    tet[j] = tetrahedra(j, i);
}

template <typename I, typename F>
void simplicial_unstructured_3d_mesh<I, F>::get_triangle(I i, I tri[]) const
{
  tri[0] = triangles(0, i);
  tri[1] = triangles(1, i);
  tri[2] = triangles(2, i);
}

template <typename I, typename F>
void simplicial_unstructured_3d_mesh<I, F>::get_edge(I i, I v[]) const
{
  v[0] = edges(0, i);
  v[1] = edges(1, i);
}

template <typename I, typename F>
void simplicial_unstructured_3d_mesh<I, F>::get_coords(I i, F coords[]) const
{
  coords[0] = vertex_coords(0, i);
  coords[1] = vertex_coords(1, i);
  coords[2] = vertex_coords(2, i);
}

template <typename I, typename F>
bool simplicial_unstructured_3d_mesh<I, F>::find_edge(const I v[2], I &i) const
{
  const auto it = edge_id_map.find( std::make_tuple(v[0], v[1]) );
  if (it == edge_id_map.end()) return false;
  else {
    i = it->second;
    assert(edges(0, i) == v[0]);
    assert(edges(1, i) == v[1]);
    return true;
  }
}

template <typename I, typename F>
bool simplicial_unstructured_3d_mesh<I, F>::find_triangle(const I v[3], I &i) const
{
  i = -1;
  const auto it = triangle_id_map.find( std::make_tuple(v[0], v[1], v[2]) );
  if (it == triangle_id_map.end()) return false;
  else {
    i = it->second;
    // fprintf(stderr, "found triangle!!! %d\n", i);
    assert(triangles(0, i) == v[0]);
    assert(triangles(1, i) == v[1]);
    assert(triangles(2, i) == v[2]);
    return true;
  }
}

template <typename I, typename F>
bool simplicial_unstructured_3d_mesh<I, F>::find_tetrahedron(const I v[4], I &i) const
{
  const auto it = tetrahedron_id_map.find( std::make_tuple(v[0], v[1], v[2], v[3]) );
  if (it == tetrahedron_id_map.end()) return false;
  else {
    i = it->second;
    assert(tetrahedra(0, i) == v[0]);
    assert(tetrahedra(1, i) == v[1]);
    assert(tetrahedra(2, i) == v[2]);
    assert(tetrahedra(3, i) == v[3]);
    return true;
  }
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_3d_mesh<I, F>::side_of(int d, I i) const
{
  if (d == 0)
    return vertex_side_of[i];
  else if (d == 1) 
    return edges_side_of[i];
  else if (d == 2)
    return triangle_side_of[i];
  else 
    return std::set<I>();
}

template <typename I, typename F>
inline void simplicial_unstructured_3d_mesh<I, F>::build_smoothing_kernel(const F sigma)
{
  const F sigma2 = sigma * sigma;
  const F limit = F(3) * sigma;

  auto neighbors = [&](I i) { return vertex_edge_vertex[i]; };
  smoothing_kernel.resize(n(0));

  for (auto i = 0; i < n(0); i ++) {
    std::set<I> set;
    const F xi[3] = {vertex_coords[i*3], vertex_coords[i*3+1], vertex_coords[i*3+2]};
    // fprintf(stderr, "i=%d, x={%f, %f}\n", i, xi[0], xi[1]);
    auto criteron = [&](I j) {
      const F xj[3] = {vertex_coords[j*3], vertex_coords[j*3+1], vertex_coords[j*3+2]};
      if (vector_dist_2norm_3(xi, xj) < limit)
        return true;
      else return false;
    };

    auto operation = [&set](I j) {set.insert(j);};
    bfs<I, std::set<I>>(i, neighbors, operation, criteron);

    auto &kernel = smoothing_kernel[i];
    for (auto k : set) {
      const F xk[3] = {vertex_coords[3*k], vertex_coords[3*k+1], vertex_coords[3*k+2]};
      const F d = vector_dist_2norm_3(xi, xk);
      const F w = std::exp(-(d*d) / (2*sigma*sigma)) / (sigma * std::sqrt(2.0 * M_PI));
      // fprintf(stderr, "d2=%f, w=%f\n", d2, w);
      kernel.push_back( std::make_tuple(k, w) );
      // fprintf(stderr, "i=%d, k=%d, %f\n", i, k, w); 
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
void simplicial_unstructured_3d_mesh<I, F>::smooth_scalar_gradient_jacobian(
    const ndarray<F>& f, const F sigma, 
    ndarray<F>& scalar, // smoothed scalar field
    ndarray<F>& grad,  // smoothed gradient field
    ndarray<F>& J) const // smoothed jacobian field
{
  const F sigma2 = sigma * sigma, 
          sigma4 = sigma2 * sigma2;

  scalar.reshape({n(0)});
  grad.reshape({3, n(0)});
  J.reshape({3, 3, n(0)});

  // for (auto i = 0; i < smoothing_kernel.size(); i ++) {
  this->parallel_for(smoothing_kernel.size(), [&](int i) {
    for (auto j = 0; j < smoothing_kernel[i].size(); j ++) {
      auto tuple = smoothing_kernel[i][j];
      const auto k = std::get<0>(tuple);
      const auto w = std::get<1>(tuple);
    
      const F d[3] = {vertex_coords[k*3] - vertex_coords[i*3], 
                      vertex_coords[k*3+1] - vertex_coords[i*3+1],
                      vertex_coords[k*3+2] - vertex_coords[i*3+2]};

      // scalar
      scalar[i] += f[k] * w;

      // gradient
      grad(0, i) += - f[k] * w * d[0] / sigma2;
      grad(1, i) += - f[k] * w * d[1] / sigma2;
      grad(2, i) += - f[k] * w * d[2] / sigma2;

      // jacobian
      J(0, 0, i) += (d[0]*d[0] / sigma2 - 1) / sigma2 * f[k] * w;
      J(0, 1, i) += d[0]*d[1] / sigma4 * f[k] * w;
      J(0, 2, i) += d[0]*d[2] / sigma4 * f[k] * w;

      J(1, 0, i) += d[0]*d[1] / sigma4 * f[k] * w;
      J(1, 1, i) += (d[1]*d[1] / sigma2 - 1) / sigma2 * f[k] * w;
      J(1, 2, i) += d[1]*d[2] / sigma4 * f[k] * w;

      J(2, 0, i) += d[2]*d[0] / sigma4 * f[k] * w;
      J(2, 1, i) += d[2]*d[1] / sigma4 * f[k] * w;
      J(2, 2, i) += (d[2]*d[2] / sigma2 - 1) / sigma2 * f[k] * w;
    }
  });
}

} // namespace ftk

///////// serialization

namespace diy {
  template <typename I, typename F> struct Serialization<ftk::simplicial_unstructured_3d_mesh<I, F>> {
    static void save(diy::BinaryBuffer& bb, const ftk::simplicial_unstructured_3d_mesh<I, F>& m) {
      diy::save(bb, m.vertex_coords);
      diy::save(bb, m.vertex_side_of);
      diy::save(bb, m.vertex_edge_vertex);
      diy::save(bb, m.edges);
      diy::save(bb, m.edges_side_of);
      diy::save(bb, m.triangles);
      diy::save(bb, m.triangle_sides);
      diy::save(bb, m.triangle_side_of);
      diy::save(bb, m.tetrahedra);
      diy::save(bb, m.tetrahedra_sides);
      diy::save(bb, m.edge_id_map);
      diy::save(bb, m.triangle_id_map);
      diy::save(bb, m.tetrahedron_id_map);
    }
   
    static void load(diy::BinaryBuffer& bb, ftk::simplicial_unstructured_3d_mesh<I, F>& m) {
      diy::load(bb, m.vertex_coords);
      diy::load(bb, m.vertex_side_of);
      diy::load(bb, m.vertex_edge_vertex);
      diy::load(bb, m.edges);
      diy::load(bb, m.edges_side_of);
      diy::load(bb, m.triangles);
      diy::load(bb, m.triangle_sides);
      diy::load(bb, m.triangle_side_of);
      diy::load(bb, m.tetrahedra);
      diy::load(bb, m.tetrahedra_sides);
      diy::load(bb, m.edge_id_map);
      diy::load(bb, m.triangle_id_map);
      diy::load(bb, m.tetrahedron_id_map);
    }
  };
} // namespace diy

#endif
