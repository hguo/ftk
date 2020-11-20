#ifndef _FTK_MESH_UNSTRUCTURED_3D_HH
#define _FTK_MESH_UNSTRUCTURED_3D_HH

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

public: // io
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_vtk_unstructured_grid() const;
  void from_vtk_unstructured_grid(vtkSmartPointer<vtkUnstructuredGrid> grid);
#endif

public: 
  void element_for(int d, std::function<void(I)> f);

public:
  void get_tetrahedron(I i, I tet[]) const;
  void get_triangle(I i, I tri[]) const;
  void get_edge(I i, I edge[]) const;

public:
  std::set<I> sides(int d, I i) const;
  std::set<I> side_of(int d, I i) const;

  void get_simplex(int d, I i, I simplex[]) const;
  void get_coords(I i, F coords[]) const;

  const ndarray<F>& get_coords() const {return vertex_coords;}

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

  initialize();
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
}

}

#endif
