#ifndef _HYPERMESH_simplicial_unstructured_extruded_2d_mesh_HH
#define _HYPERMESH_simplicial_unstructured_extruded_2d_mesh_HH

#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_unstructured_extruded_2d_mesh : public object { // extruded from 
  simplicial_unstructured_extruded_2d_mesh(const simplicial_unstructured_2d_mesh<I, F>& m_) : m(m_) {extrude();}

  size_t n(int d) const;
  size_t n_ordinal(int d) const { return m.n(d); }
  size_t n_interval(int d) const;

  bool is_ordinal(int d, I i) const { const auto type = mod(i, n(d)); return type < m.n(d); }

  I flat_vertex_id(I i) const { return mod(i, m.n(0)); }
  I flat_vertex_time(I i) const { return i / m.n(0); }
  I extruded_vertex_id(I i, bool t=true) { return t ? i + m.n(0) : i; }

  int face_type(I i) const;

public: // element iteration
  void element_for(int d, std::function<void(I)> f, 
      int nthreads=std::thread::hardware_concurrency()) const;

  void element_for_ordinal(int d, int t, std::function<void(I)> f,
      int nthreads=std::thread::hardware_concurrency()) const;
  
  void element_for_interval(int d, int t, std::function<void(I)> f,
      int nthreads=std::thread::hardware_concurrency()) const;

public: // mesh access
  std::set<I> sides(int d, I i) const;
  std::set<I> side_of(int d, I i) const;

  void get_simplex(int d, I i, I verts[]) const;
  void get_coords(I i, F coords[]) const;
  
protected:
  void extrude();

  static I mod(I val, I m);

private:
  const simplicial_unstructured_2d_mesh<I, F>& m;

  ftk::ndarray<I> tets, // {4, n(3)}
                  tris, // {3, n(2)}
                  edges;// {2, n(1)}

  ftk::ndarray<I> tets_sides, // {4, n(3)}
                  tris_sides; // {3, n(2)}
};

/////
template <typename I, typename F>
I simplicial_unstructured_extruded_2d_mesh<I, F>::mod(I v, I m)
{
  I mod = v % m;
  if (v < 0) mod += m;
  return mod;
}

template <typename I, typename F>
void simplicial_unstructured_extruded_2d_mesh<I, F>::extrude() 
{
  tets.reshape({4, n(3)});
  
  // tets
  for (auto i = 0; i < m.n(2); i ++) {
    I tri[3];
    m.get_triangle(i, tri);
    std::sort(tri, tri+3);

    // tets
    tets(0, i) = tri[0]; // type I: bottom tet
    tets(1, i) = tri[1];
    tets(2, i) = tri[2];
    tets(3, i) = tri[2] + m.n(0);

    tets(0, i+m.n(2)) = tri[0]; // type II: center tet
    tets(1, i+m.n(2)) = tri[1];
    tets(2, i+m.n(2)) = tri[1] + m.n(0);
    tets(3, i+m.n(2)) = tri[2] + m.n(0);

    tets(0, i+2*m.n(2)) = tri[0]; // type III: top tet
    tets(1, i+2*m.n(2)) = tri[0] + m.n(0);
    tets(2, i+2*m.n(2)) = tri[1] + m.n(0);
    tets(3, i+2*m.n(2)) = tri[2] + m.n(0);
  }

  // unique triangles
  tris.reshape({3, n(2)}); 
  for (auto i = 0; i < m.n(2); i ++) { // per prism
    I tri[3];
    m.get_triangle(i, tri);
    
    tris(0, i) = tri[0]; // type I: base triangle
    tris(1, i) = tri[1];
    tris(2, i) = tri[2];

    std::sort(tri, tri+3); // the following extrusions need sorted vertices
    tris(0, i + m.n(2)) = tri[0]; // type II: prism lower triangle
    tris(1, i + m.n(2)) = tri[1];
    tris(2, i + m.n(2)) = tri[2] + m.n(0);
    
    tris(0, i + m.n(2)*2) = tri[0]; // type III: prism upper triangle
    tris(1, i + m.n(2)*2) = tri[1] + m.n(0);
    tris(2, i + m.n(2)*2) = tri[2] + m.n(0);
  }
  for (auto i = 0; i < m.n(1); i ++) { // per edge 
    I v[2];
    m.get_edge(i, v);
    if (v[0] > v[1]) std::swap(v[0], v[1]);

    tris(0, i + m.n(2)*3) = v[0]; // type IV: edge lower triangle
    tris(1, i + m.n(2)*3) = v[1];
    tris(2, i + m.n(2)*3) = v[1] + m.n(0);
    
    tris(0, i + m.n(2)*3 + m.n(1)) = v[0]; // type V: edge upper triangle
    tris(1, i + m.n(2)*3 + m.n(1)) = v[0] + m.n(0);
    tris(2, i + m.n(2)*3 + m.n(1)) = v[1] + m.n(0);
  }

  // unique edges
  edges.reshape({2, n(1)});
  for (auto i = 0; i < m.n(1); i ++) {
    I v[2];
    m.get_edge(i, v);
    if (v[0] > v[1]) std::swap(v[0], v[1]);

    edges(0, i) = v[0];
    edges(1, i) = v[1];
    
    edges(0, i+m.n(1)) = v[0];
    edges(1, i+m.n(1)) = v[1] + m.n(0);
    
    edges(0, i+m.n(1)*2) = v[0];
    edges(1, i+m.n(1)*2) = v[0] + m.n(0);
  }
}

template <typename I, typename F>
size_t simplicial_unstructured_extruded_2d_mesh<I, F>::n(int d) const
{
  if (d == 0) {
    return m.n(0); // unique vertices in each interval
  } else if (d == 1) {
    return m.n(1) * 3; // unique edges in each interval
  } else if (d == 2) {
    return m.n(2) // the bottom triangle of the prism
      + m.n(2) * 2 // two triangles inside prisms
      + m.n(1) * 2;// unique triangles extruded from each edge
  } else if (d == 3) { 
    return 3 * m.n(2); // 3 tets per prism
  } else return 0;
}

template <typename I, typename F>
size_t simplicial_unstructured_extruded_2d_mesh<I, F>::n_interval(int d) const
{
  if (d == 1) return m.n(1)*2;
  else if (d == 2) return m.n(2)*2 + m.n(1)*2;
  else if (d == 3) return m.n(2)*3;
  else return 0;
}

template <typename I, typename F>
void simplicial_unstructured_extruded_2d_mesh<I, F>::element_for(int d, std::function<void(I)> f, int nthreads) const
{
  parallel_for(n(d), nthreads, [&](int i) {f(i);});
  // for (auto i = 0; i < n(d); i ++)
  //   f(i);
}

template <typename I, typename F>
void simplicial_unstructured_extruded_2d_mesh<I, F>::element_for_ordinal(int d, int t, std::function<void(I)> f, int nthreads) const
{
  parallel_for(n_ordinal(d), nthreads, [&](int i) {
      f(i + t * n(d));});

  // for (auto i = 0; i < n_ordinal(d); i ++)
  //   f(i + t * n(d));
}

template <typename I, typename F>
void simplicial_unstructured_extruded_2d_mesh<I, F>::element_for_interval(int d, int t, std::function<void(I)> f, int nthreads) const
{
  parallel_for(n_interval(d), nthreads, [&](int i) {
    const auto id = i + n_ordinal(d) + t * n(d);
    f(id);
  });

  // for (auto i = 0; i < n_interval(d); i ++) {
  //   const auto id = i + n_ordinal(d) + t * n(d);
  //   f(id);
    // if (d == 2) side_of(2, id);
  // }
}

template <typename I, typename F>
void simplicial_unstructured_extruded_2d_mesh<I, F>::get_simplex(int d, I k, I verts[]) const
{
  const I i = mod(k, n(d)), t = std::floor(double(k) / n(d));
  const I offset = t * n(0);
  // fprintf(stderr, "d=%d, k=%d, i=%d, t=%d, offset=%d\n", d, k, i, t, offset);
  if (d == 0) {
    verts[0] = i + offset;
  } else if (d == 1) {
    verts[0] = edges(0, i) + offset;
    verts[1] = edges(1, i) + offset;
  } else if (d == 2) {
    verts[0] = tris(0, i) + offset;
    verts[1] = tris(1, i) + offset;
    verts[2] = tris(2, i) + offset;
  } else if (d == 3) {
    verts[0] = tets(0, i) + offset;
    verts[1] = tets(1, i) + offset;
    verts[2] = tets(2, i) + offset;
    verts[3] = tets(3, i) + offset;
  }
}

template <typename I, typename F>
void simplicial_unstructured_extruded_2d_mesh<I, F>::get_coords(I i, F coords[]) const
{
  const I k = flat_vertex_id(i),
          t = flat_vertex_time(i);
  m.get_coords(k, coords);
  coords[2] = t;
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_extruded_2d_mesh<I, F>::sides(int d, I k) const 
{
  const I i = mod(k, n(d)), t = std::floor((double)k / n(d));
  std::set<I> results;

  I v[4];
  get_simplex(3, i, v);

  if (d == 3) {
    const int type = i / m.n(2);
    // fprintf(stderr, "tet_type=%d\n", type);
    if (type == 0) { // bottom 0-1-2-2'
      // 0-1-2, type I
      // 0-1-2', type II
      // 1-2-2', type IV
      // 0-2-2', type IV
      {
        const I ot[3] = {v[0], v[1], v[2]}; // triangle in the orignal mesh
        I otid; // triangle id in the original mesh
        bool found = m.find_triangle(ot, otid);
        assert(found);

        results.insert( otid + t*n(2) );
        results.insert( otid + t*n(2) + m.n(2) );
      }
      {
        const I oe[2] = {v[1], v[2]};
        I oeid; // edge id in the orignal mesh
        bool found = m.find_edge(oe, oeid);
        assert(found);
        results.insert( oeid + m.n(2)*3 + t*n(2));
      }
      {
        const I oe[2] = {v[0], v[2]};
        I oeid; // edge id in the orignal mesh
        bool found = m.find_edge(oe, oeid);
        assert(found);
        results.insert( oeid + m.n(2)*3 + t*n(2));
      }
    } else if (type == 1) { // center 0-1-1'-2'
      // 0-1-2', type II
      // 0-1'-2', type III
      // 0-1-1', type IV
      // 1-1'-2', type V
      {
        const I ot[3] = {v[0], v[1], mod(v[3], m.n(0))}; // triangle in the orignal mesh
        I otid; // triangle id in the original mesh
        bool found = m.find_triangle(ot, otid);
        assert(found);

        results.insert( otid + m.n(2) + t*n(2) );
        results.insert( otid + m.n(2)*2 + t*n(2) );
      }
      {
        const I oe[2] = {v[0], v[1]};
        I oeid; // edge id in the orignal mesh
        bool found = m.find_edge(oe, oeid);
        assert(found);
        results.insert( oeid + m.n(2)*3 + t*n(2) );
      }
      {
        const I oe[2] = {v[1], mod(v[3], n(0))}; 
        I oeid; // edge id in the orignal mesh
        bool found = m.find_edge(oe, oeid);
        assert(found);
        results.insert( oeid + m.n(2)*3 + m.n(1) + t*n(2) );
      }
    } else if (type == 2) { // top 0-0'-1'-2'
      // 0-1-2, t+1, type I
      // 0-1'-2', type III
      // 0-0'-1', type V
      // 0-0'-2', type V
      {
        const I ot[3] = {v[0], mod(v[2], m.n(0)), mod(v[3], m.n(0))}; // triangle in the orignal mesh
        I otid; // triangle id in the original mesh
        bool found = m.find_triangle(ot, otid);
        assert(found);

        results.insert( otid + (t+1)*n(2) );
        results.insert( otid + m.n(2)*2 + t*n(2) );
      }
      {
        const I oe[2] = {v[0], mod(v[2], n(0))}; 
        I oeid; // edge id in the orignal mesh
        bool found = m.find_edge(oe, oeid);
        assert(found);
        results.insert( oeid + m.n(2)*3 + m.n(1) + t*n(2) );
      }
      {
        const I oe[2] = {v[0], mod(v[3], m.n(0))}; 
        I oeid; // edge id in the orignal mesh
        bool found = m.find_edge(oe, oeid);
        assert(found);
        results.insert( oeid + m.n(2)*3 + m.n(1) + t*n(2) );
      }
    }
  }

  return results;
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_extruded_2d_mesh<I, F>::side_of(int d, I k) const 
{
  const I i = mod(k, n(d)), t = std::floor((double)k / n(d));
  // fprintf(stderr, "k=%d, i=%d, t=%d, n(d)=%d\n", k, i, t, n(d));
  std::set<I> results;

  if (d == 2) {
    int v[3];
    get_simplex(2, k, v);

    // fprintf(stderr, "k=%d, i=%d, t=%d, v=%d, %d, %d\n", k, i, t, v[0], v[1], v[2]);
    // fprintf(stderr, "face=%d,%d type=%d\n", k, i, face_type(k));

    // determine triangle type by id
    const I v0[3] = {mod(v[0], m.n(0)), mod(v[1], m.n(0)), mod(v[2], m.n(0))};
    if (i < 3*m.n(2)) {
      I otid; // triangle id in the original mesh
      bool found = m.find_triangle(v0, otid); 
      assert(found);
      
      if (i < m.n(2)) { // type I: bottom tet
        // t-1: 0-0'-1'-2' (type III tet)
        // t:   0-1-2-2' (type I tet)
        results.insert( otid + 2*m.n(2) + (t-1)*n(3) );
        results.insert( otid            + t    *n(3) );
      } else if (i < 2*m.n(2)) { // type II: prism lower triangle
        // t:   0-1-2-2' (type I tet)
        // t:   0-1-1'-2' (type II tet)
        results.insert( otid          + t*n(3) );
        results.insert( otid + m.n(2) + t*n(3) );
      } else { // if (i < 3*m.n(2)) { // type III: prism upper triangle
        // t:   0-1-1'-2'  (type II tet)
        // t:   0-0'-1'-2' (type III tet)
        results.insert( otid +   m.n(2) + t*n(3) );
        results.insert( otid + 2*m.n(2) + t*n(3) );
      }
    } else {
      I oeid; // edge id in the original mesh
      I ve[2] = {v0[0], v0[2]};
      bool found = m.find_edge(ve, oeid);
      assert(found);

      const bool lower = i < 3*m.n(2) + m.n(1);

      const auto triangles = m.side_of(1, oeid);
      for (const auto triangle : triangles) {
        I vc[3];
        m.get_triangle(triangle, vc);

        I pos = -1;
        for (int j = 0; j < 3; j ++)
          if (vc[j] < ve[0]) pos = 0; // small
          else if (vc[j] > ve[0] && vc[j] < ve[1]) pos = 1; // middle
          else if (vc[j] > ve[1]) pos = 2; // large
        assert(pos != -1);

        if (lower) { //type IV: edge lower triangle
          if (pos == 0) 
            results.insert(triangle + t*n(3)); // bottom tet
          else if (pos == 1)
            results.insert(triangle + t*n(3)); // bottom tet
          else if (pos == 2)
            results.insert(triangle + m.n(2) + t*n(3)); // center tet
        } else { // type V: edge upper triangle
          if (pos == 0) 
            results.insert(triangle + m.n(2) + t*n(3)); // center tet
          else if (pos == 1)
            results.insert(triangle + 2*m.n(2) + t*n(3)); // top tet
          else if (pos == 2)
            results.insert(triangle + 2*m.n(2) + t*n(3)); // top tet
        }
      }

#if 0
      int tids[2] = {-1}; // triangle ids in array
      int j = 0;
      for (auto tri : triangles) {
        tids[j++] = tri;

        I ve1[2];
        m.get_edge(oeid, ve1);

        int vv[3];
        m.get_triangle(tri, vv);
        fprintf(stderr, "ve=%d, %d, oeid=%d, ve1=%d, %d, vv=%d, %d, %d, tri=%d\n", 
            ve[0], ve[1], oeid, 
            ve1[0], ve1[1], 
            vv[0], vv[1], vv[2], tri);
      }
#endif

#if 0
      // find triangle v0-v1-v2, denote the other triangle as vn-v0-v1
      if (i < 3*m.n(2) + m.n(1)) { // type IV: edge lower triangle
        // type II tet for the large triangle
        // type I tet for the small triangle
        if (large_triangle >= 0) 
          results.insert( large_triangle + m.n(2) + t*n(3) );
        if (small_triangle >= 0) 
          results.insert( small_triangle          + t*n(3) );
      } else { // type V: edge upper triangle
        // type III tet for the large triangle
        // type II tet for the small triangle
        if (large_triangle >= 0)
          results.insert( large_triangle + 2*m.n(2) + t*n(3) );
        if (small_triangle >= 0) 
          results.insert( small_triangle +   m.n(2) + t*n(3) );
      }
#endif
    }


    // for debug only
    bool fail = false;
    for (const auto tet : results) {
      const auto s = sides(3, tet);
      if (s.find(k) == s.end()) fail = true;
    }

    if (fail) {
      auto tm = [&](I i, int d=0) { return I(std::floor(double(i)/n(d))); };

      fprintf(stderr, "n(d)=%zu, %zu, %zu, %zu\n", n(0), n(1), n(2), n(3));
      fprintf(stderr, "face=%d(%d): %d(%d), %d(%d), %d(%d)\n", 
          k, tm(k, 2), v[0], tm(v[0]), v[1], tm(v[1]), v[2], tm(v[2]));
      for (const auto tet : results) {
        int v[4];
        get_simplex(3, tet, v);
        fprintf(stderr, "--tet=%d(%d): %d(%d), %d(%d), %d(%d), %d(%d)\n",
          tet, tm(tet, 3), v[0], tm(v[0]), v[1], tm(v[1]), v[2], tm(v[2]), v[3], tm(v[3]));
        fprintf(stderr, "----sides:");
        for (auto s : sides(3, tet))
          fprintf(stderr, "%d ", s);
        fprintf(stderr, "\n");
      }
      assert(false);
    }
  }

  return results;
}

template <typename I, typename F>
int simplicial_unstructured_extruded_2d_mesh<I, F>::face_type(I k) const
{
  const I i = mod(k, n(2));

  if (i < m.n(2)) return 0;  // I
  else if (i < 2*m.n(2)) return 1; // II
  else if (i < 3*m.n(2)) return 2; // III 
  else if (i < 3*m.n(2) + m.n(1)) return 3; // IV
  else if (i < 3*m.n(2) + 2*m.n(1)) return 4; // V
  else return -1;
}

}

#endif
