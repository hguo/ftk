#ifndef _HYPERMESH_SIMPLEX_2D_HH
#define _HYPERMESH_SIMPLEX_2D_HH

#include <ftk/ftk_config.hh>
#include <ftk/hypermesh/ndarray.hh>
#include <iostream>
#include <vector>

namespace hypermesh {

struct simplex_2d_mesh;

struct simplex_2d_mesh_element {

};

struct simplex_2d_mesh { // 2D triangular mesh
  friend class simplex_2d_mesh_element;
  typedef simplex_2d_mesh_element iterator;

  simplex_2d_mesh() {}

  simplex_2d_mesh(
      const std::vector<double>& coords, // coordinates of vertices; the dimension of the array is 2 * n_vertices
      const std::vector<size_t>& triangles); // vertex id of triangles; the dimension of the array is 3 * n_triangles

  simplex_2d_mesh(
      const ndarray<double>& coords_, // 2 * n_vertices
      const ndarray<size_t>& triangles_) // 3 * n_triangles
    : vertex_coords(coords_), triangles(triangles_) {}

  // dimensionality of the mesh
  int nd() const {return 2;}

  // numer of d-dimensional elements
  size_t n(int d) const;

  void build_edges();

private:
  ndarray<double> vertex_coords; // 2 * n_vertices
  std::vector<std::vector<size_t>> vertex_side_of;

  ndarray<size_t> triangles; // 3 * n_triangles

  ndarray<size_t> edges; // 2 * n_edges
  ndarray<size_t> edges_side_of; // 2 * n_edges
};

/////////

simplex_2d_mesh::simplex_2d_mesh(const std::vector<double> &coords_, const std::vector<size_t> &triangles_)
{
  vertex_coords.copy_vector(coords_);
  vertex_coords.reshape({2, coords_.size()/2});
  
  triangles.copy_vector(triangles_);
  triangles.reshape({3, triangles_.size()/3});
}

size_t simplex_2d_mesh::n(int d) const
{
  if (d == 0) return vertex_coords.dim(1);
  else if (d == 1) { // TODO
    return 0;
  } else if (d == 2)
    return triangles.dim(2);
}

void simplex_2d_mesh::build_edges()
{
  typedef std::tuple<size_t, size_t> edge_t;
  std::set<edge_t> unique_edges;

  std::auto convert_edge = [](edge_t e) {
    if (std::get<0>(e) > std::get<1>(e))
      return std::make_tuple(std::get<1>(e), std::get<0>(e));
    else return e;
  }

  for (auto i = 0; i < triangles.dim(1); i ++) {
    unique_edges.insert(std::make_tuple(triangles(0, i), triangles(1, i)));
    unique_edges.insert(std::make_tuple(triangles(1, i), triangles(2, i)));
    unique_edges.insert(std::make_tuple(triangles(2, i), triangles(3, i)));
  }

  edges.reshape(2, unique_edges.size());

  size_t i = 0;
  for (const auto e : unique_edges) {
    edges(0, i) = std::get<0>(e);
    edges(1, i) = std::get<1>(e);
    i ++;
  }
}

}

#endif
