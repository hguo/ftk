#ifndef _FTK_UNSTRUCTURED_MESH_2D_HH
#define _FTK_UNSTRUCTURED_MESH_2D_HH

#include <ftk/config.hh>
#include <ftk/algorithms/bfs.hh>
#include <ftk/mesh/simplicial_unstructured_mesh.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/sign_det.hh>
#include <ftk/utils/bcast.hh>

#if FTK_HAVE_METIS
#include <metis.h>
#endif

namespace ftk {

template <typename I, typename F> struct point_locator_2d;

template <typename I=int>
struct hash_compare_tuple_i2 {
  hash_compare_tuple_i2() {}
  hash_compare_tuple_i2(const hash_compare_tuple_i2&) {}
  ~hash_compare_tuple_i2() {}

  bool equal(const std::tuple<I, I>& a, const std::tuple<I, I>& b) const { return a == b; }
  size_t hash(const std::tuple<I, I>& x) const {
    I a = std::get<0>(x), b = std::get<1>(x);
    return a * b * 2654435761;
  };
};

template <typename I=int>
struct hash_compare_tuple_i3 {
  hash_compare_tuple_i3() {}
  hash_compare_tuple_i3(const hash_compare_tuple_i3&) {}
  ~hash_compare_tuple_i3() {}

  bool equal(const std::tuple<I, I, I>& a, const std::tuple<I, I, I>& b) const { return a == b; }
  size_t hash(const std::tuple<I, I, I>& x) const {
    I a = std::get<0>(x), b = std::get<1>(x), c = std::get<2>(x);
    return a * b * c * 2654435761;
  };
};

template <typename I=int, typename F=double>
struct simplicial_unstructured_2d_mesh : // 2D triangular mesh
  public simplicial_unstructured_mesh<I, F>
{ 
  friend class point_locator_2d<I, F>;
  friend class diy::Serialization<simplicial_unstructured_2d_mesh<I, F>>;

  simplicial_unstructured_2d_mesh() {}

  simplicial_unstructured_2d_mesh(
      const std::vector<F>& coords, // coordinates of vertices; the dimension of the array is 2 * n_vertices
      const std::vector<I>& triangles); // vertex id of triangles; the dimension of the array is 3 * n_triangles

  simplicial_unstructured_2d_mesh(
      const ndarray<F>& coords_, // 2 * n_vertices
      const ndarray<I>& triangles_) // 3 * n_triangles
    : vertex_coords(coords_), triangles(triangles_) {build_triangles(); build_edges(); build_partition();}

  // dimensionality of the mesh
  int nd() const {return 2;}

  // numer of d-dimensional elements
  size_t n(int d, bool part=false) const;

  void build_edges();
  void build_triangles();

  void build_smoothing_kernel_cached(F sigma);
  std::string default_smoothing_kernel_filename(F sigma) const;

  bool has_smoothing_kernel() const { return smoothing_kernel.size() > 0; }
  void build_smoothing_kernel(F sigma);
  const std::vector<std::vector<std::tuple<I, F>>> &get_smoothing_kernel() const { return smoothing_kernel; }
  F get_smoothing_kernel_size() const { return sigma; }
  void smooth_scalar_gradient_jacobian(
      const ndarray<F>& f, 
      // const F sigma,
      ndarray<F>& fs, // smoothed scalar field
      ndarray<F>& g,  // smoothed gradient field
      ndarray<F>& j   // smoothed jacobian field
  ) const; 
  ndarray<F> smooth_scalar(const ndarray<F>& f) const;

  ndarray<F> scalar_gradient(const ndarray<F>& f) const;
  ndarray<F> vector_gradient(const ndarray<F>& f) const;
  void scalar_gradient_jacobian(const ndarray<F>& f, ndarray<F>& g, ndarray<F>& j) const;

  template <typename T> bool eval(const ndarray<T>& f, const F x[], T val[]) const; // will not work for 2.5d meshes

public: // distributed parallel
  void build_partition();
  const std::vector<I>& get_part_triangles() const { return part_triangles; }

  std::shared_ptr<simplicial_unstructured_2d_mesh<I, F>> generate_submesh() const;

public: // io
  void from_vtu(const std::string filename);
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkUnstructuredGrid> to_vtu() const;
  void from_vtu(vtkSmartPointer<vtkUnstructuredGrid> grid);
#endif
  void to_vtu(const std::string& filename) const;

  template <typename T>
  void array_to_vtu(const std::string& filename, const std::string& varname, const ndarray<T>&) const;
  
  template <typename T>
  void array_to_vtu(const std::string& filename, 
      const std::vector<std::string>& varnames, 
      const std::vector<ndarray<T>>& arrays) const;

  void write_smoothing_kernel(const std::string& filename);
  bool read_smoothing_kernel(const std::string& filename);

public: // element iteration
  void element_for(int d, std::function<void(I)> f);

public: // mesh access
  std::set<I> sides(int d, I i) const;
  std::set<I> side_of(int d, I i) const;

  // std::set<I> side_of2(const I v[2]) const;

  I get_simplex(int d, I i, I v[], bool part=false) const; // note that the output verts are always in gid; always return gid
  I get_triangle(I i, I tri[], bool part=false) const; // always return gid
  I get_edge(I i, I edge[], bool part=false) const; // always return gid
  void get_coords(I i, F coords[], bool part=false) const;

  bool find_simplex(int d, const I v[], I& i) const;
  bool find_edge(const I v[2], I& i) const;
  bool find_triangle(const I v[3], I& i) const;

  const ndarray<F>& get_coords() const {return vertex_coords;}
  const ndarray<I>& get_triangles() const {return triangles;}
  const ndarray<I>& get_edges() const {return edges;}
  const std::vector<std::set<I>>& get_vertex_triangles() const {return vertex_triangles;}
  
  const std::set<I>& get_vertex_edge_vertex(I i) const {return vertex_edge_vertex[i];}
  
  int get_triangle_chi(I i) const { return triangles_chi[i] ? -1 : 1; }

  // virtual int ncoords() const { return vertex_coords.dim(0); } // { return 2; }
  // virtual int ncoords() const { return 3; }

public: // point locator and misc
  I nearest(F x[]) const; // locate which point is nearest to x
  // I locate(F x[]) const; // locate which triangle contains x
  I locate(F x[], F mu[]) const { return locator->locate(x, mu); } // locate which triangle contains x and get the barycentric coordinates of x

  std::shared_ptr<point_locator_2d<I, F>> get_locator() const { return locator; }

protected:
  mutable std::shared_ptr<point_locator_2d<I, F>> locator;

protected: // mesh connectivities
  ndarray<F> vertex_coords; // 2 * n_vertices
  std::vector<std::set<I>> vertex_side_of;
  std::vector<std::set<I>> vertex_edge_vertex;
  std::vector<std::set<I>> vertex_triangles;

  ndarray<I> edges; // 2 * n_edges
  ndarray<I> edges_side_of; // 2 * n_edges
  // std::map<std::tuple<I, I>, std::set<I>> edges_side_of;

  ndarray<I> triangles; // 3 * n_triangles
  ndarray<I> triangle_edges; // 3 * n_triangles
  std::vector<bool> triangles_chi; // chiralities of triangles, 1 means +1; 0 means -1
  
  std::vector<std::set<I>> triangle_edge_triangles;

public: // additional mesh info
#if FTK_HAVE_TBB
  tbb::concurrent_hash_map<std::tuple<I, I>, int, hash_compare_tuple_i2<I>> edge_id_map;
  tbb::concurrent_hash_map<std::tuple<I, I, I>, int, hash_compare_tuple_i3<I>> triangle_id_map;
#else
  std::map<std::tuple<I, I>, int> edge_id_map;
  std::map<std::tuple<I, I, I>, int> triangle_id_map;
#endif

public: // partition
  bool is_partial() const { return partial; }
  I lid2gid(int d, I) const; // translate ID in the local partition to global ID
  I gid2lid(int d, I) const; // translate global ID to local ID

protected: // parallel partition
  bool partial = false;
  std::vector<I> part_triangles, part_edges, part_vertices;
  std::map<I/*gid*/, I/*lid*/> part_triangles_gid, part_edges_gid, part_vertices_gid;

private:
  F sigma; // smoothing kernel size
  std::vector<std::vector<std::tuple<I/*vert*/, F/*weight*/>>> smoothing_kernel;
};


/////////

template <typename I, typename F>
template <typename T>
bool simplicial_unstructured_2d_mesh<I, F>::eval(const ndarray<T>& f, const F x[], T val[]) const
{
  F mu[3];
  I tid = this->locator->locate(x, mu);
  if (tid < 0) {
    return false;
  } else {
    I tri[3];
    this->get_simplex(2, tid, tri);
   
    if (f.multicomponents() == 0) { // scalar field
      T fs[3];
      for (int i = 0; i < 3; i ++)
        fs[i] = f[tri[i]];
      val[0] = lerp_s2(fs, mu);
      return true;
    } else if (f.multicomponents() == 1) { // vector field
      const int d = f.dim(0);
      T fs[3][d];

      for (int i = 0; i < 3; i ++)
        for (int j = 0; j < d; j ++)
          fs[i][j] = f(j, tri[i]);

      for (int j = 0; j < d; j ++)
        val[j] = lerp_s2(fs[j], mu);
      return true;
    } else // not supported
      return false;
  }
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::build_smoothing_kernel_cached(F sigma)
{
  const auto f = default_smoothing_kernel_filename(sigma);
  if (file_exists(f))
    read_smoothing_kernel(f);
  else {
    build_smoothing_kernel(sigma);
    write_smoothing_kernel(f);
  }
}

template <typename I, typename F>
std::string simplicial_unstructured_2d_mesh<I, F>::default_smoothing_kernel_filename(F sigma) const
{
  unsigned int h0 = triangles.hash();
  unsigned int h1 = murmurhash2(&sigma, sizeof(F), h0);

  std::stringstream ss;
  ss << "ftk.smoothing.kernels.2d." 
     << std::hex << h1;
  return ss.str();
}

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
size_t simplicial_unstructured_2d_mesh<I, F>::n(int d, bool part) const
{
  if (part) {
    if (d == 0) return part_vertices.size(); 
    else if (d == 1) return part_edges.size(); 
    else if (d == 2) return part_triangles.size();
    else return 0;
  } else {
    if (d == 0) return vertex_coords.dim(1);
    else if (d == 1) {
      return edges.dim(1);
    } else if (d == 2)
      return triangles.dim(1);
    else return 0;
  }
}

template <typename I, typename F>
I simplicial_unstructured_2d_mesh<I, F>::get_simplex(int d, I i, I v[], bool part) const
{
  if (d == 1) 
    return get_edge(i, v, part);
  else if (d == 2) 
    return get_triangle(i, v, part);
  else 
    return -1; // assert false;
}

template <typename I, typename F>
bool simplicial_unstructured_2d_mesh<I, F>::find_simplex(int d, const I v[3], I& i) const
{
  if (d == 1) return find_edge(v, i);
  else if (d == 2) return find_triangle(v, i);
  return false;
}

template <typename I, typename F>
bool simplicial_unstructured_2d_mesh<I, F>::find_edge(const I v_[2], I &i) const
{
  I v[2] = {v_[0], v_[1]};
  if (v[0] > v[1]) std::swap(v[0], v[1]);

#if FTK_HAVE_TBB
  typename tbb::concurrent_hash_map<std::tuple<I, I>, int, hash_compare_tuple_i2<I>>::const_accessor it;
  bool found = edge_id_map.find(it, std::make_tuple(v[0], v[1]));
  if (found) {
    i = it->second;
    assert(edges(0, i) == v[0]);
    assert(edges(1, i) == v[1]);
    return true;
  } else 
    return false;
#else
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
#endif
}

template <typename I, typename F>
bool simplicial_unstructured_2d_mesh<I, F>::find_triangle(const I v_[3], I &i) const
{
  int v[3] = {v_[0], v_[1], v_[2]};
  std::sort(v, v+3);

#if FTK_HAVE_TBB
  typename tbb::concurrent_hash_map<std::tuple<I, I, I>, int, hash_compare_tuple_i3<I>>::const_accessor it;
  bool found = triangle_id_map.find(it, std::make_tuple(v[0], v[1], v[2]));
  if (found) {
    i = it->second;
    assert(triangles(0, i) == v[0]);
    assert(triangles(1, i) == v[1]);
    assert(triangles(2, i) == v[2]);
    return true;
  } else 
    return false;
#else
  const auto it = triangle_id_map.find( std::make_tuple(v[0], v[1], v[2]) );
  if (it == triangle_id_map.end()) return false;
  else {
    i = it->second;
    assert(triangles(0, i) == v[0]);
    assert(triangles(1, i) == v[1]);
    assert(triangles(2, i) == v[2]);
    return true;
  }
#endif
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::build_triangles()
{
  vertex_triangles.resize(n(0));
  triangles_chi.resize(n(2), 0);

  for (auto i = 0; i < n(2); i ++) {
    I v[3], order[3];
    get_triangle(i, v);
    // fprintf(stderr, "id=%d, v=%d, %d, %d\n", 
    //     i, v[0], v[1], v[2]);
    const int nswaps = nswaps_bubble_sort<3, I>(v, order);
    // std::sort(v, v+3);
    
    if (nswaps % 2 == 0) triangles_chi[i] = 0;
    else triangles_chi[i] = 1;
    // fprintf(stderr, "id=%d, nswaps=%d\n", i, nswaps);

    for (auto j = 0; j < 3; j ++)
      triangles(j, i) = v[j];

    // triangle_id_map[ std::make_tuple(v[0], v[1], v[2]) ] = i;
    triangle_id_map.insert( { std::make_tuple(v[0], v[1], v[2]), i } );

    for (int j = 0; j < 3; j ++)
      vertex_triangles[v[j]].insert(i);
  }
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::build_edges()
{
  // std::map<std::tuple<I, I>, std::set<I>> map_edges_side_of;
  std::map<I, std::set<I>> map_edges_side_of;

  int edge_count = 0;
  auto add_edge = [&](I tid, I k, I v0, I v1) {
    if (v0 > v1) std::swap(v0, v1);
    const auto edge = std::make_tuple(v0, v1);
    I id;
#if FTK_HAVE_TBB
    typename tbb::concurrent_hash_map<std::tuple<I, I>, int, hash_compare_tuple_i2<I>>::const_accessor it;
    bool found = edge_id_map.find(it, edge);
    if (!found) {
      id = edge_count ++;
      edge_id_map.insert( {edge, id} );
    } else 
      id = it->second;
#else
    if (edge_id_map.find(edge) == edge_id_map.end()) {
      id = edge_count ++;
      edge_id_map[edge] = id;
    } else 
      id = edge_id_map[edge];
#endif

    map_edges_side_of[id].insert(tid);
    triangle_edges(k, tid) = id;

    // fprintf(stderr, "adding edge #%d (%d, %d) from triangle %d\n", id, v0, v1, tid);
    // if (v0 == 0 && v1 == 1) {
    //   fprintf(stderr, "{0, 1}, tid=%d\n", tid);
    //   exit(1);
    // }
  };

  // fprintf(stderr, "triangles_88078=%d, %d, %d\n", triangles(0, 88078), triangles(1, 88078), triangles(2, 88078)); 
  triangle_edges.reshape(3, n(2));
  for (auto i = 0; i < n(2); i ++) {
    add_edge(i, 0, triangles(0, i), triangles(1, i));
    add_edge(i, 1, triangles(1, i), triangles(2, i));
    add_edge(i, 2, triangles(2, i), triangles(0, i));
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

  // triangle_edge_triangles
  triangle_edge_triangles.resize(n(2));
  for (int i = 0; i < n(2); i ++) {
    const I t0 = edges_side_of(0, i), 
            t1 = edges_side_of(1, i);
    if (t0 >= 0 && t1 >= 0) {
      triangle_edge_triangles[t0].insert(t1);
      triangle_edge_triangles[t1].insert(t0);
    }
  }
}

template <typename I, typename F>
inline void simplicial_unstructured_2d_mesh<I, F>::build_smoothing_kernel(const F sigma)
{
  fprintf(stderr, "building smoothing kernel...\n");
  this->sigma = sigma;

  const F sigma2 = sigma * sigma;
  const F limit = F(3) * sigma;

  auto neighbors = [&](I i) { return vertex_edge_vertex[i]; };
  smoothing_kernel.resize(n(0));

  // for (auto i = 0; i < n(0); i ++) {
  object::parallel_for(n(0), [&](int i) {
    std::set<I> set;
    // const F xi[2] = {vertex_coords[i*2], vertex_coords[i*2+1]};
    const F xi[2] = {vertex_coords(0, i), vertex_coords(1, i)};

    auto criteron = [&](I j) {
      // const F xj[2] = {vertex_coords[j*2], vertex_coords[j*2+1]};
      const F xj[2] = {vertex_coords(0, j), vertex_coords(1, j)};
      if (vector_dist_2norm_2(xi, xj) < limit) {
        // fprintf(stderr, "%d (%f, %f) add %d (%f, %f)\n", i, xi[0], xi[1], j, xj[0], xj[1]);
        // fprintf(stderr, "i=%d, x={%f, %f}\n", i, xi[0], xi[1]);
        return true;
      }
      else return false;
    };

    auto operation = [&set](I j) {set.insert(j);};
    // fprintf(stderr, "bfs for node %d...\n", i);
    bfs<I, std::set<I>>(i, neighbors, operation, criteron);
    // fprintf(stderr, "bfs for node %d done.\n", i);

    auto &kernel = smoothing_kernel[i];
    for (auto k : set) {
      // const F xk[2] = {vertex_coords[2*k], vertex_coords[2*k+1]};
      const F xk[2] = {vertex_coords(0, k), vertex_coords(1, k)};
      const F d = vector_dist_2norm_2(xi, xk);
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
  });
}

template <typename I, typename F>
ndarray<F> simplicial_unstructured_2d_mesh<I, F>::smooth_scalar(const ndarray<F>& f) const
{
  const F sigma2 = sigma * sigma, 
          sigma4 = sigma2 * sigma2;

  ndarray<F> scalar;
  scalar.reshape({n(0)});

  // for (auto i = 0; i < smoothing_kernel.size(); i ++) {
  this->parallel_for(smoothing_kernel.size(), [&](int i) {
    for (auto j = 0; j < smoothing_kernel[i].size(); j ++) {
      auto tuple = smoothing_kernel[i][j];
      const auto k = std::get<0>(tuple);
      const auto w = std::get<1>(tuple);

      scalar[i] += f[k] * w;
    }
  });

  return scalar;
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::smooth_scalar_gradient_jacobian(
    const ndarray<F>& f, // const F sigma, 
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
  this->parallel_for(smoothing_kernel.size(), [&](int i) {
    for (auto j = 0; j < smoothing_kernel[i].size(); j ++) {
      auto tuple = smoothing_kernel[i][j];
      const auto k = std::get<0>(tuple);
      const auto w = std::get<1>(tuple);
    
      // const F d[2] = {vertex_coords[k*2] - vertex_coords[i*2], 
      //                 vertex_coords[k*2+1] - vertex_coords[i*2+1]};
      const F d[2] = {vertex_coords(0, k) - vertex_coords(0, i),
                      vertex_coords(1, k) - vertex_coords(1, i)};
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
ndarray<F> simplicial_unstructured_2d_mesh<I, F>::vector_gradient(const ndarray<F>& vector) const
{
  // compute cellwise gradient
  ndarray<F> cellgrad;
  cellgrad.reshape({2, vector.dim(0), n(2)});
  cellgrad.set_multicomponents(2);
  for (int i = 0; i < n(2); i ++) {
    I tri[3];
    get_triangle(i, tri);
   
    F X[3][2]; 
    for (int j = 0; j < 3; j ++) {
      const I k = tri[j];
      X[j][0] = vertex_coords(0, k);
      X[j][1] = vertex_coords(1, k);
    }
    
    for (int c = 0; c < vector.dim(0); c ++) {
      F f[3], gradf[2];
      for (int j = 0; j < 3; j ++) {
        f[j] = vector(c, tri[j]);
      }

      gradient_2dsimplex2(X, f, gradf);
      cellgrad(0, c, i) = gradf[0];
      cellgrad(1, c, i) = gradf[1];
    }
  }
  
  // compute pointwise gradient
  ndarray<F> grad;
  grad.reshape({2, vector.dim(0), n(0)});
  grad.set_multicomponents(2);
  for (int i = 0; i < n(0); i ++) {
    const auto tris = vertex_triangles[i];
    const auto ntris = tris.size();

    for (const auto tri : tris) {
      for (int c = 0; c < vector.dim(0); c ++) {
        grad(0, c, i) += cellgrad(0, c, tri);
        grad(1, c, i) += cellgrad(1, c, tri);
      }
    }

    for (int c = 0; c < vector.dim(0); c ++) {
      grad(0, c, i) = grad(0, c, i) / ntris;
      grad(1, c, i) = grad(0, c, i) / ntris;
    }
  }
  
  return grad;
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::scalar_gradient_jacobian(const ndarray<F>& f, ndarray<F>& g, ndarray<F>& j) const
{
  g = scalar_gradient(f);
  j = vector_gradient(g);
}

template <typename I, typename F>
ndarray<F> simplicial_unstructured_2d_mesh<I, F>::scalar_gradient(const ndarray<F>& scalar) const
{
  // compute cellwise gradient
  ndarray<F> cellgrad;
  cellgrad.reshape({2, n(2)});
  cellgrad.set_multicomponents();
  for (int i = 0; i < n(2); i ++) {
    I tri[3];
    get_triangle(i, tri);
    
    F X[3][2], f[3], gradf[2];
    for (int j = 0; j < 3; j ++) {
      const I k = tri[j];
      X[j][0] = vertex_coords(0, k);
      X[j][1] = vertex_coords(1, k);
      f[j] = scalar[k];
    }

    gradient_2dsimplex2(X, f, gradf);
    cellgrad(0, i) = gradf[0];
    cellgrad(1, i) = gradf[1];
  }
  
  // compute pointwise gradient
  ndarray<F> grad;
  grad.reshape({2, n(0)});
  grad.set_multicomponents();
  for (int i = 0; i < n(0); i ++) {
    F gradf[2] = {0};
    const auto tris = vertex_triangles[i];
    const auto ntris = tris.size();
    for (const auto tri : tris) {
      gradf[0] += cellgrad(0, tri);
      gradf[1] += cellgrad(1, tri);
    }

    grad(0, i) = gradf[0] / ntris;
    grad(1, i) = gradf[1] / ntris;
  }
  
  return grad;
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_2d_mesh<I, F>::sides(int d, I i) const
{
  std::set<I> results;
  if (d == 1) {
    results.insert( edges(0, i) );
    results.insert( edges(1, i) );
  } else if (d == 2) {
    results.insert( triangle_edges(0, i) );
    results.insert( triangle_edges(1, i) );
    results.insert( triangle_edges(2, i) );
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

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::from_vtu(const std::string filename)
{
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkXMLUnstructuredGridReader::New();
  reader->SetFileName(filename.c_str());
  reader->Update();
  vtkSmartPointer<vtkUnstructuredGrid> grid = reader->GetOutput();
  from_vtu(grid);
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

#if FTK_HAVE_VTK
template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::from_vtu(vtkSmartPointer<vtkUnstructuredGrid> grid)
{
  if (this->is_root_proc()) {
    vtkIdType ncells = grid->GetNumberOfCells();
    std::vector<int> m_triangles;
    for (vtkIdType i = 0; i < ncells; i ++) {
      vtkSmartPointer<vtkCell> cell = grid->GetCell(i);
      if (cell->GetCellType() == VTK_TRIANGLE) {
        vtkIdType v[3] = {cell->GetPointId(0), cell->GetPointId(1), cell->GetPointId(2)};
        // std::sort(v, v+3);
        for (int j = 0; j < 3; j ++)
          m_triangles.push_back(v[j]);
      }
    }
    triangles.reshape({3, m_triangles.size()/3});
    triangles.from_vector(m_triangles);

    vtkIdType npts = grid->GetNumberOfPoints();
    vertex_coords.reshape({3, size_t(npts)});
    for (vtkIdType i = 0; i < npts; i ++) {
      double x[3] = {0};
      grid->GetPoint(i, x);
      vertex_coords(0, i) = x[0];
      vertex_coords(1, i) = x[1];
      vertex_coords(2, i) = x[2];
    }
  }

  diy::mpi::bcastv(this->comm, triangles, this->get_root_proc());
  diy::mpi::bcastv(this->comm, vertex_coords, this->get_root_proc());

  build_triangles();
  build_edges();
  build_partition();
}

template <typename I, typename F>
vtkSmartPointer<vtkUnstructuredGrid> simplicial_unstructured_2d_mesh<I, F>::to_vtu() const
{
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> pts = vtkPoints::New();
  pts->SetNumberOfPoints(n(0));

  for (int i=0; i<n(0); i++) {
    F coords[3] = {0};
    get_coords(i, coords);
    pts->SetPoint(i, coords[0], coords[1], coords[2]);
    // pts->SetPoint(i, vertex_coords[i*2], vertex_coords[i*2+1], 0); 
  }

  for (int i=0; i<n(2); i ++) {
    vtkIdType ids[3] = {triangles[i*3], triangles[i*3+1], triangles[i*3+2]};
    grid->InsertNextCell(VTK_TRIANGLE, 3, ids);
  }

  grid->SetPoints(pts);
  return grid;
}
#endif

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::to_vtu(const std::string& filename) const
{
#if FTK_HAVE_VTK
  auto grid = to_vtu();
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInputData( grid );
  writer->Write();
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

template <typename I, typename F>
template <typename T>
void simplicial_unstructured_2d_mesh<I, F>::array_to_vtu(
    const std::string& filename, 
    const std::vector<std::string>& varnames, 
    const std::vector<ndarray<T>>& arrays) const
{
#if FTK_HAVE_VTK
  auto grid = to_vtu();

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

template <typename I, typename F>
template <typename T>
void simplicial_unstructured_2d_mesh<I, F>::array_to_vtu(const std::string& filename, const std::string& varname, const ndarray<T>& array) const
{
#if FTK_HAVE_VTK
  auto grid = to_vtu();
  auto data = array.to_vtk_data_array(varname);

  grid->GetPointData()->AddArray(data);
  grid->GetPointData()->SetActiveScalars(varname.c_str());

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInputData( grid );
  writer->Write();
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::write_smoothing_kernel(const std::string& f)
{
  fprintf(stderr, "writing smoothing kernel to %s\n", f.c_str());
  diy::serializeToFile(sigma, smoothing_kernel, f);
}

template <typename I, typename F>
bool simplicial_unstructured_2d_mesh<I, F>::read_smoothing_kernel(const std::string& f)
{
  fprintf(stderr, "reading smoothing kernel from %s\n", f.c_str());
  bool succ = diy::unserializeFromFile(f, sigma, smoothing_kernel);
  return succ;
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::element_for(int d, std::function<void(I)> f)
{
  for (auto i = 0; i < n(d); i ++)
    f(i);
}

template <typename I, typename F>
I simplicial_unstructured_2d_mesh<I, F>::get_triangle(I i, I tri[], bool part) const
{
  if (part) 
    return get_triangle(lid2gid(2, i), tri, false);
  else {
    tri[0] = triangles[i*3];
    tri[1] = triangles[i*3+1];
    tri[2] = triangles[i*3+2];
    return i;
  }
#if 0
  tri[0] = triangles(0, i);
  tri[1] = triangles(1, i);
  tri[2] = triangles(2, i);
#endif
}

template <typename I, typename F>
I simplicial_unstructured_2d_mesh<I, F>::get_edge(I i, I v[], bool part) const
{
  if (part)
    return get_edge(lid2gid(1, i), v, false);
  else {
    v[0] = edges(i*2);
    v[1] = edges(i*2+1);
    return i;
  }
  // v[0] = edges(0, i);
  // v[1] = edges(1, i);
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::get_coords(I i, F coords[], bool part) const
{
  if (part)
    get_coords(lid2gid(0, i), coords, false);
  else {
    coords[0] = vertex_coords(0, i);
    coords[1] = vertex_coords(1, i);
    if (vertex_coords.dim(0) > 2)
      coords[2] = vertex_coords(2, i);
    else 
      coords[2] = F(0);
  }
}

template <typename I, typename F>
I simplicial_unstructured_2d_mesh<I, F>::nearest(F x[]) const
{
  // brute-force search, for now
  F mindist = std::numeric_limits<F>::max();
  I imindist;
  for (int i = 0; i < n(0); i ++) {
    // F y[2] = {vertex_coords[2*i], vertex_coords[2*i+1]};
    F y[2] = {vertex_coords(0, i), vertex_coords(1, i)};
    F dist = vector_dist_2norm_2(x, y);
    if (mindist > dist) { 
      mindist = dist;
      imindist = i;
    }
  }
  return imindist;
}

template <typename I, typename F>
I simplicial_unstructured_2d_mesh<I, F>::lid2gid(int d, I lid) const 
{
  if (d == 2) { // triangles
    return part_triangles[lid];
  } else if (d == 1) { // edges
    return part_edges[lid];
  } else if (d == 0) { // vertices
    return part_vertices[lid];
  }
  return -1;
}

template <typename I, typename F>
I simplicial_unstructured_2d_mesh<I, F>::gid2lid(int d, I gid) const 
{
  if (d == 2) { // triangles
    return part_triangles_gid.at(gid);
  } else if (d == 1) { // edges
    return part_edges_gid.at(gid);
  } else if (d == 0) { // vertices
    return part_vertices_gid.at(gid);
  }
  return -1;
}

template <typename I, typename F>
void simplicial_unstructured_2d_mesh<I, F>::build_partition()
{
  // TODO: use parmetis
#if FTK_HAVE_METIS
  if (this->comm.size() == 1) return;
    
  idx_t nverts = n(2), ncon = 1, nparts = this->comm.size();
  std::vector<idx_t> xadj, adj, part; // (nverts, 0);
  idx_t objval;

  if (this->comm.rank() == 0) {
    fprintf(stderr, "building partitions...\n");

    for (int i = 0; i < n(2); i ++) {
      xadj.push_back(adj.size());
      for (const auto t : triangle_edge_triangles[i])
        adj.push_back(t);
    }
    xadj.push_back(adj.size());
    part.resize(nverts, 0);

    int rtn = METIS_PartGraphKway(
        &nverts, 
        &ncon, 
        &xadj[0],
        &adj[0],
        NULL, 
        NULL, 
        NULL, 
        &nparts,
        NULL,
        NULL,
        NULL,
        &objval,
        &part[0]);

    fprintf(stderr, "partition succ.\n");
  }
  diy::mpi::bcastv(this->comm, part);

  // part_triangles
  part_triangles.clear();
  part_triangles_gid.clear();
  std::set<I> set_part_edges, set_part_verts;
  
  for (I i = 0; i < n(2); i ++) {
    if (part[i] == this->comm.rank()) {
      // fprintf(stderr, "rank=%d, tri=%d\n", this->comm.rank(), i);
      const I lid = part_triangles.size();
      part_triangles.push_back(i);
      part_triangles_gid[lid] = i;

      for (const auto eid : sides(2, i))
        set_part_edges.insert(eid);

      I tri[3];
      get_triangle(i, tri);
      for (int k = 0; k < 3; k ++)
        set_part_verts.insert(tri[k]);
    }
  }

  // fprintf(stderr, "stage1\n");

  int n_local_edges = 0;
  for (const I eid : set_part_edges) {
    const I local_eid = n_local_edges ++;
    part_edges.push_back(eid);
    part_edges_gid[local_eid] = eid;
  }

  int n_local_verts = 0;
  for (const I vid : set_part_edges) {
    const I local_vid = n_local_verts ++;
    part_vertices.push_back(vid);
    part_vertices_gid[local_vid] = vid;
  }

  partial = true;

  fprintf(stderr, "partition: rank=%d, #tri=%zu, #edge=%zu, #verts=%zu\n", this->comm.rank(), 
      part_triangles.size(), 
      part_edges.size(), 
      part_vertices.size());
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_METIS);
#endif
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

} // namespace ftk

///////// serialization

namespace diy {
  template <typename I, typename F> struct Serialization<ftk::simplicial_unstructured_2d_mesh<I, F>> {
    static void save(diy::BinaryBuffer& bb, const ftk::simplicial_unstructured_2d_mesh<I, F>& m) {
      diy::save(bb, m.vertex_coords);
      diy::save(bb, m.vertex_side_of);
      diy::save(bb, m.vertex_edge_vertex);
      diy::save(bb, m.edges);
      diy::save(bb, m.edges_side_of);
      diy::save(bb, m.triangles);
      diy::save(bb, m.triangle_edges);
      diy::save(bb, m.edge_id_map);
      diy::save(bb, m.triangle_id_map);
    }
   
    static void load(diy::BinaryBuffer& bb, ftk::simplicial_unstructured_2d_mesh<I, F>& m) {
      diy::load(bb, m.vertex_coords);
      diy::load(bb, m.vertex_side_of);
      diy::load(bb, m.vertex_edge_vertex);
      diy::load(bb, m.edges);
      diy::load(bb, m.edges_side_of);
      diy::load(bb, m.triangles);
      diy::load(bb, m.triangle_edges);
      diy::load(bb, m.edge_id_map);
      diy::load(bb, m.triangle_id_map);
    }
  };
} // namespace diy


#endif
