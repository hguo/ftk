#ifndef _HYPERMESH_SIMPLEX_2D_HH
#define _HYPERMESH_SIMPLEX_2D_HH

#include <ftk/ftk_config.hh>
#include <ftk/ndarray.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/algorithms/bfs.hh>
#include <set>
#include <iostream>
#include <vector>

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

  void build_edges();

  void build_smoothing_kernel(F sigma);
  void smooth_scalar_field(const ndarray<F> &f);

  std::set<I> sides(int d, I i);
  std::set<I> side_of(int d, I i);

private: // mesh connectivities
  ndarray<F> vertex_coords; // 2 * n_vertices
  std::vector<std::set<I>> vertex_side_of;

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
  else if (d == 1) { // TODO FIXME
    return 0;
  } else if (d == 2)
    return triangles.dim(2);
  else return 0;
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
  for (const auto e : unique_edges) {
    auto v0 = edges(0, i) = std::get<0>(e);
    auto v1 = edges(1, i) = std::get<1>(e);
    vertex_side_of[v0].insert(i);
    vertex_side_of[v1].insert(i);
    i ++;
  }
}

template <typename I, typename F>
inline void simplex_2d_mesh<I, F>::build_smoothing_kernel(const F sigma)
{
  const F limit2 = std::pow(F(3) * sigma, F(2));

  auto neighbors = [&](I i) {
    std::set<I> results;
    for (auto edge : side_of(0, i))
      for (auto side : sides(1, edge))
        results.insert(side);
    return results;
  };

  smoothing_kernel.resize(n(0));

  for (auto i = 0; i < n(0); i ++) {
    std::set<I> set;
    const F xi[2] = {vertex_coords(0, i), vertex_coords(1, i)};
    // fprintf(stderr, "i=%d, x={%f, %f}\n", i, xi[0], xi[1]);
    auto criteron = [&](I j) {
      const F xj[2] = {vertex_coords(0, j), vertex_coords(1, j)};
      if (vector_dist_2norm_2(xi, xj) < limit2)
        return true;
      else return false;
    };

    auto operation = [&set](I j) {set.insert(j);};
    bfs<I, std::set<I>>(i, neighbors, operation, criteron);

    auto &kernel = smoothing_kernel[i];
    for (auto k : set) {
      const F xk[2] = {vertex_coords(0, k), vertex_coords(1, k)};
      const F d2 = vector_dist_2norm_2(xi, xk);
      const F w = std::exp(-d2 / (sigma*sigma)) / (sigma * std::sqrt(2.0 * M_PI));
      // fprintf(stderr, "d2=%f, w=%f\n", d2, w);
      kernel.push_back( std::make_tuple(k, w) );
    }

    // normalization
    F sum = 0;
    for (int k = 0; k < kernel.size(); k ++)
      sum += std::get<1>(kernel[k]);
    for (int k = 0; k < kernel.size(); k ++) {
      std::get<1>(kernel[k]) /= sum;
      // fprintf(stderr, "i=%d, k=%d, %f\n", i, k, std::get<1>(kernel[k]));// kernel.size());
    }
  }
}

template <typename I, typename F>
inline void simplex_2d_mesh<I, F>::smooth_scalar_field(const ndarray<F>& f)
{
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

}

#endif
