#ifndef _HYPERMESH_SIMPLEX_2D_EXTRUSION_MESH_HH
#define _HYPERMESH_SIMPLEX_2D_EXTRUSION_MESH_HH

#include <ftk/mesh/simplex_2d_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplex_2d_extrusion_mesh { // extruded from 
  simplex_2d_extrusion_mesh(const simplex_2d_mesh<I, F>& m_) : m(m_) {extrude();}

  size_t n(int d) const;

  I flat_vertex_id(I i) { return i % m.n(0); }
  I extrude_vertex_id(I i, bool t=true) { return t ? i + m.n(0) : i; }

public: // element iteration
  void element_for(int d, std::function<void(I)> f);

public: // mesh access
  std::set<I> sides(int d, I i);
  std::set<I> side_of(int d, I i);

  void get_simplex(I i, I verts[]) const;
  
protected:
  void extrude();

private:
  const simplex_2d_mesh<I, F>& m;

  ftk::ndarray<I> tets, // {4, n(3)}
                  tris, // {3, n(2)}
                  edges;// {2, n(1)}

  ftk::ndarray<T> tets_sides, // {4, n(3)}
                  tris_sides; // {3, n(2)}
};

/////
template <typename I, typename F>
void simplex_2d_extrusion_mesh<I, F>::extrude() 
{
  tets.reshape({4, n(3)});
  tris.reshape({3, n(2)});
  edges.reshape({2, n(1)});

  for (auto i = 0; i < m.n(2); i ++) {
    I tri[3];
    m.get_triangle(i, tri);
    std::sort(tri, tri+3);

    // tets
    tets(0, i*3) = tri[0];
    tets(1, i*3) = tri[1];
    tets(2, i*3) = tri[2] + m.n(0);
    tets(3, i*3) = tri[1] + m.n(0);
    
    tets(0, i*3+1) = tri[0];
    tets(1, i*3+1) = tri[1];
    tets(2, i*3+1) = tri[2];
    tets(3, i*3+1) = tri[2] + m.n(0);

    tets(0, i*3+2) = tri[0];
    tets(1, i*3+2) = tri[1] + m.n(0);
    tets(2, i*3+2) = tri[2] + m.n(0);
    tets(3, i*3+2) = tri[0] + m.n(0);

    // unique triangles
    tris(0, i*9) = tri[0];
    tris(1, i*9) = tri[1];
    tris(2, i*9) = tri[2];
    
    tris(0, i*9+1) = tri[0];
    tris(1, i*9+1) = tri[1];
    tris(2, i*9+1) = tri[1] + m.n(0);

    tris(0, i*9+2) = tri[0];
    tris(1, i*9+2) = tri[1] + m.n(0);
    tris(2, i*9+2) = tri[0] + m.n(0);

    tris(0, i*9+3) = tri[1];
    tris(1, i*9+3) = tri[2];
    tris(2, i*9+3) = tri[2] + m.n(0);

    tris(0, i*9+4) = tri[1];
    tris(1, i*9+4) = tri[2] + m.n(0);
    tris(2, i*9+4) = tri[1] + m.n(0);
    
    tris(0, i*9+5) = tri[0];
    tris(1, i*9+5) = tri[2];
    tris(2, i*9+5) = tri[2] + m.n(0);

    tris(0, i*9+6) = tri[0];
    tris(1, i*9+6) = tri[2] + m.n(0);
    tris(2, i*9+6) = tri[0] + m.n(0);
    
    tris(0, i*9+7) = tri[0];
    tris(1, i*9+7) = tri[1];
    tris(2, i*9+7) = tri[2] + m.n(0);
    
    tris(0, i*9+8) = tri[0];
    tris(1, i*9+8) = tri[1] + m.n(0);
    tris(2, i*9+8) = tri[2] + m.n(0);
    
    // unique edges
    edges(0, i*9) = tri[0];
    edges(1, i*9) = tri[1];
    
    edges(0, i*9+1) = tri[1];
    edges(1, i*9+1) = tri[2];
    
    edges(0, i*9+2) = tri[2];
    edges(1, i*9+2) = tri[0];
    
    edges(0, i*9+3) = tri[0];
    edges(1, i*9+3) = tri[0] + m.n(0);
    
    edges(0, i*9+4) = tri[1];
    edges(1, i*9+4) = tri[1] + m.n(0);
    
    edges(0, i*9+5) = tri[2];
    edges(1, i*9+5) = tri[2] + m.n(0);
    
    edges(0, i*9+6) = tri[0];
    edges(1, i*9+6) = tri[1] + m.n(0);
    
    edges(0, i*9+7) = tri[1];
    edges(1, i*9+7) = tri[2] + m.n(0);
    
    edges(0, i*9+8) = tri[0];
    edges(1, i*9+8) = tri[2] + m.n(0);
  }
}

template <typename I, typename F>
size_t simplex_2d_extrusion_mesh<I, F>::n(int d) const
{
  if (d == 0) {
    return m.n(0); // unique vertices in each interval
  } else if (d == 1) { 
    return m.n(2) * 9; // unique edges in each interval
  } else if (d == 2) { 
    return m.n(2) * 9; // unique triangles in each interval
  } else if (d == 3) { 
    return 3 * m.n(2); // 3 tets per prism
  } else return 0;
}

}

#endif
