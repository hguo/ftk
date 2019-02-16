#ifndef _FTK_CURVE2TUBE_HH
#define _FTK_CURVE2TUBE_HH

#include <tuple>
#include <vector>
#include <ftk/numeric/vector_dot_product.hh>
#include <ftk/numeric/vector_assignment.hh>
#include <ftk/numeric/vector_normalization.hh>

namespace ftk {

template <typename R> // float or double
std::tuple<std::vector<unsigned int>/*indices*/, std::vector<R>/*vertices*/, std::vector<R>/*normals*/, std::vector<R>/*colors*/>
curve2tube(const std::vector<R>& curve_vertices, const std::vector<R>& curve_colors, int npatches, R radius)
{
  std::vector<unsigned int> tube_indices;
  std::vector<R> tube_vertices, tube_normals, tube_colors;
  const int curve_nvertices = curve_vertices.size() / 3;

  if (curve_nvertices < 2) return std::make_tuple(tube_indices, tube_vertices, tube_normals, tube_colors);

  R N0[3] = {0};
  for (int j = 1; j < curve_nvertices; j ++) {
    R color[3] = {curve_colors[j*3], curve_colors[j*3+1], curve_colors[j*3+2]};
    R P0[3] = {curve_vertices[(j-1)*3], curve_vertices[(j-1)*3+1], curve_vertices[(j-1)*3+2]}, 
      P[3]  = {curve_vertices[j*3], curve_vertices[j*3+1], curve_vertices[j*3+2]};

    R T[3];
    vector_assign_subtraction3(T, P, P0);
    vector_normalization2_3(T);

    R N[3] = {-T[1], T[0], R(0)};
    vector_normalization2_3(N);
    if (!std::isnormal(vector_2norm_3(N)))
      vector_assign3(N, R(1), R(0), R(0));

    R B[3];
    cross_product(N, T, B);

    if (j > 1) {
      R n0 = ftk::vector_dot_product3(N0, N);
      R b0 = ftk::vector_dot_product3(N0, B);
      
      R N1[3];
      vector_assign_weighted_addition3(N1, n0, N, b0, B);
      vector_normalization2_3(N1);
      vector_assign3(N, N1);

      cross_product(N, T, B);
      vector_normalization2_3(B);
    }
    ftk::vector_assign3(N0, N);

    const int nIteration = (j==1) ? 2 : 1; 
    for (int k = 0; k < nIteration; k ++) {
      for (int p = 0; p < npatches; p ++) {
        R angle = p * 2 * M_PI / npatches;

        R normal[3];
        vector_assign_weighted_addition3(normal, cos(angle), N, sin(angle), B);
        vector_normalization2_3(normal);

        R offset[3];
        vector_assign_scalar_multiplication3(offset, radius, normal);

        R coords[3];
        if (k == 0 && j == 1) vector_assign_addition3(coords, P0, offset);
        else vector_assign_addition3(coords, P, offset);

        tube_vertices.push_back(coords[0]);
        tube_vertices.push_back(coords[1]);
        tube_vertices.push_back(coords[2]);

        tube_normals.push_back(normal[0]);
        tube_normals.push_back(normal[1]);
        tube_normals.push_back(normal[2]);

        tube_colors.push_back(color[0]); 
        tube_colors.push_back(color[1]); 
        tube_colors.push_back(color[2]);
      }
    }

    for (int p=0; p<npatches; p++) {
      const int n = tube_vertices.size()/3; 
      const int pn = (p + 1) % npatches; 
      tube_indices.push_back(n-npatches+p); 
      tube_indices.push_back(n-npatches-npatches+pn); 
      tube_indices.push_back(n-npatches-npatches+p); 
      tube_indices.push_back(n-npatches+p); 
      tube_indices.push_back(n-npatches+pn); 
      tube_indices.push_back(n-npatches-npatches+pn); 
    }
  }

  return std::make_tuple(tube_indices, tube_vertices, tube_normals, tube_colors);
}

}

#endif
