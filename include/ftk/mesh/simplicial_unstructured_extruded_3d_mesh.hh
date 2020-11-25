#ifndef _HYPERMESH_simplicial_unstructured_extruded_3d_mesh_HH
#define _HYPERMESH_simplicial_unstructured_extruded_3d_mesh_HH

#include <ftk/mesh/simplicial_unstructured_3d_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_unstructured_extruded_3d_mesh : public object { // extruded from 
  simplicial_unstructured_extruded_3d_mesh(const simplicial_unstructured_3d_mesh<I, F>& m_) : m(m_) {extrude();}

  size_t n(int d) const;
  size_t n_ordinal(int d) const { return m.n(d); }
  size_t n_interval(int d) const;

  bool is_ordinal(int d, I i) const {
    I k = mod(i, n(3));
    return k < n_ordinal(d);
  }

  I flat_vertex_id(I i) const { return mod(i, m.n(0)); }
  I flat_vertex_time(I i) const { return i / m.n(0); }
  I extruded_vertex_id(I i, bool t=true) { return t ? i + m.n(0) : i; }

  int face_type(I i) const;

public: // element iteration
  void element_for(int d, std::function<void(I)> f, 
      int xl = FTK_XL_NONE, int nthreads=std::thread::hardware_concurrency(), bool affinity = false) const;

  void element_for_ordinal(int d, int t, std::function<void(I)> f,
      int xl = FTK_XL_NONE, int nthreads=std::thread::hardware_concurrency(), bool affinity = false) const;
  
  void element_for_interval(int d, int t, std::function<void(I)> f,
      int xl = FTK_XL_NONE, int nthreads=std::thread::hardware_concurrency(), bool affinity = false) const;

public: // mesh access
  std::set<I> sides(int d, I i) const;
  std::set<I> side_of(int d, I i) const;

  void get_simplex(int d, I i, I verts[]) const;
  void get_coords(I i, F coords[]) const;
  
protected:
  void extrude();

  static I mod(I val, I m);

private:
  const simplicial_unstructured_3d_mesh<I, F>& m;

  ftk::ndarray<I> pents,// {5, n(4)}
                  tets, // {4, n(3)}
                  tris, // {3, n(2)}
                  edges;// {2, n(1)}

  ftk::ndarray<I> pents_sides,// {5, n(4)}
                  tets_sides, // {4, n(3)}
                  tris_sides; // {3, n(2)}
};

/////
template <typename I, typename F>
I simplicial_unstructured_extruded_3d_mesh<I, F>::mod(I v, I m)
{
  I mod = v % m;
  if (v < 0) mod += m;
  return mod;
}

template <typename I, typename F>
void simplicial_unstructured_extruded_3d_mesh<I, F>::extrude() 
{
  // pents
  pents.reshape({5, n(4)});
  for (auto i = 0; i < m.n(3); i ++) {
    I tet[4];
    m.get_tetrahedron(i, tet);

    pents(0, i) = tet[0]; // type I:   0 1 2 3 3'
    pents(1, i) = tet[1];
    pents(2, i) = tet[2];
    pents(3, i) = tet[3];
    pents(4, i) = tet[3] + m.n(0);

    pents(0, i+m.n(3)) = tet[0]; // type II:  0 1 2 2'3'
    pents(1, i+m.n(3)) = tet[1];
    pents(2, i+m.n(3)) = tet[2];
    pents(3, i+m.n(3)) = tet[2] + m.n(0);
    pents(4, i+m.n(3)) = tet[3] + m.n(0);

    pents(0, i+2*m.n(3)) = tet[0]; // type III: 0 1 1'2'3'
    pents(1, i+2*m.n(3)) = tet[1];
    pents(2, i+2*m.n(3)) = tet[1] + m.n(0);
    pents(3, i+2*m.n(3)) = tet[2] + m.n(0);
    pents(4, i+2*m.n(3)) = tet[3] + m.n(0);

    pents(0, i+3*m.n(3)) = tet[0]; // type IV:  0 0'1'2'3'
    pents(1, i+3*m.n(3)) = tet[0] + m.n(0);
    pents(2, i+3*m.n(3)) = tet[1] + m.n(0);
    pents(3, i+3*m.n(3)) = tet[2] + m.n(0);
    pents(4, i+3*m.n(3)) = tet[3] + m.n(0);
  }

  // unique tets
  tets.reshape({4, n(3)});
  for (auto i = 0; i < m.n(3); i ++) { // per tetrahedra of the base
    I tet[4];
    m.get_tetrahedron(i, tet);

    tets(0, i) = tet[0]; // 0 1 2 3, "base"
    tets(1, i) = tet[1];
    tets(2, i) = tet[2];
    tets(3, i) = tet[3];

    tets(0, i+m.n(3)) = tet[0]; // 0 1 2 3'
    tets(1, i+m.n(3)) = tet[1];
    tets(2, i+m.n(3)) = tet[2];
    tets(3, i+m.n(3)) = tet[3] + m.n(0);

    tets(0, i+2*m.n(3)) = tet[0]; // 0 1 2'3'
    tets(1, i+2*m.n(3)) = tet[1];
    tets(2, i+2*m.n(3)) = tet[2] + m.n(0);
    tets(3, i+2*m.n(3)) = tet[3] + m.n(0);

    tets(0, i+3*m.n(3)) = tet[0]; // 0 1'2'3'
    tets(1, i+3*m.n(3)) = tet[1] + m.n(0);
    tets(2, i+3*m.n(3)) = tet[2] + m.n(0);
    tets(3, i+3*m.n(3)) = tet[3] + m.n(0);
  }
  for (auto i = 0; i < m.n(2); i ++) { // per triangular side of the base
    I tri[3];
    m.get_triangle(i, tri);

    tets(0, i+4*m.n(3)) = tri[0]; // 0 1 2 2'
    tets(1, i+4*m.n(3)) = tri[1];
    tets(2, i+4*m.n(3)) = tri[2];
    tets(3, i+4*m.n(3)) = tri[2] + m.n(0);

    tets(0, i+4*m.n(3)+m.n(2)) = tri[0]; // 0 1 1'2'
    tets(1, i+4*m.n(3)+m.n(2)) = tri[1];
    tets(2, i+4*m.n(3)+m.n(2)) = tri[1] + m.n(0);
    tets(3, i+4*m.n(3)+m.n(2)) = tri[2] + m.n(0);
    
    tets(0, i+4*m.n(3)+2*m.n(2)) = tri[0]; // 0 0'1'2'
    tets(1, i+4*m.n(3)+2*m.n(2)) = tri[0] + m.n(0);
    tets(2, i+4*m.n(3)+2*m.n(2)) = tri[1] + m.n(0);
    tets(3, i+4*m.n(3)+2*m.n(2)) = tri[2] + m.n(0);
  }

  // unique triangles
  tris.reshape({3, n(2)});
  for (auto i = 0; i < m.n(2); i ++) { // per triangular side of the base
    I tri[3];
    m.get_triangle(i, tri);

    tris(0, i) = tri[0]; // 0 1 2, the base triangle
    tris(1, i) = tri[1];
    tris(2, i) = tri[2];

    tris(0, i+m.n(2)) = tri[0]; // 0 1 2'
    tris(1, i+m.n(2)) = tri[1];
    tris(2, i+m.n(2)) = tri[2] + m.n(0);
    
    tris(0, i+2*m.n(2)) = tri[0]; // 0 1'2'
    tris(1, i+2*m.n(2)) = tri[1] + m.n(0);
    tris(2, i+2*m.n(2)) = tri[2] + m.n(0);
  }
  for (auto i = 0; i < m.n(1); i ++) { // per edge of the base tet
    I edge[2];
    m.get_edge(i, edge);

    tris(0, i+3*m.n(2)) = edge[0]; // 0 1 1'
    tris(1, i+3*m.n(2)) = edge[1];
    tris(2, i+3*m.n(2)) = edge[1] + m.n(0);

    tris(0, i+3*m.n(2)+m.n(1)) = edge[0]; // 0 0'1'
    tris(1, i+3*m.n(2)+m.n(1)) = edge[0] + m.n(0);
    tris(2, i+3*m.n(2)+m.n(1)) = edge[1] + m.n(0);
  }

  // unique edges
  edges.reshape({2, n(1)});
  for (auto i = 0; i < m.n(1); i ++) {
    I v[2];
    m.get_edge(i, v);

    edges(0, i) = v[0]; // base edge
    edges(1, i) = v[1];

    edges(0, i+m.n(1)) = v[0]; // 0 1'
    edges(1, i+m.n(1)) = v[1] + m.n(0);
  }
  for (auto i = 0; i < m.n(0); i ++) {
    edges(0, i+2*m.n(1)) = i;
    edges(1, i+2*m.n(1)) = i + m.n(0);
  }

  fprintf(stderr, "4d mesh initialized, #pent=%zu, #tet=%zu, #tri=%zu, #edge=%zu, #vert=%zu\n", 
      n(4), n(3), n(2), n(1), n(0));
}

template <typename I, typename F>
size_t simplicial_unstructured_extruded_3d_mesh<I, F>::n(int d) const
{
  if (d == 0) {
    return m.n(0); // unique vertices in each interval
  } else if (d == 1) {
    return 2 * m.n(1) + m.n(0);
  } else if (d == 2) {
    return 3 * m.n(2) + 2 * m.n(1);
  } else if (d == 3) { 
    return 4 * m.n(3) + 3 * m.n(2);
  } else if (d == 4) {
    return 4 * m.n(3); // 4 pents per prism
  } else return 0;
}

template <typename I, typename F>
size_t simplicial_unstructured_extruded_3d_mesh<I, F>::n_interval(int d) const
{
  if (d == 1) return m.n(1) + m.n(0);
  else if (d == 2) return 2 * m.n(2) + 2 * m.n(1);
  else if (d == 3) return 3 * m.n(3) + 3 * m.n(2);
  else if (d == 4) return 4 * m.n(3);
  else return 0;
}

template <typename I, typename F>
void simplicial_unstructured_extruded_3d_mesh<I, F>::element_for(int d, std::function<void(I)> f, 
    int xl, int nthreads, bool affinity) const
{
  parallel_for(n(d), [&](int i) {f(i);}, xl, nthreads, affinity);
  // for (auto i = 0; i < n(d); i ++)
  //   f(i);
}

template <typename I, typename F>
void simplicial_unstructured_extruded_3d_mesh<I, F>::element_for_ordinal(int d, int t, std::function<void(I)> f, 
    int xl, int nthreads, bool affinity) const
{
  parallel_for(n_ordinal(d), [&](int i) {
      f(i + t * n(d));}, xl, nthreads, affinity);

  // for (auto i = 0; i < n_ordinal(d); i ++)
  //   f(i + t * n(d));
}

template <typename I, typename F>
void simplicial_unstructured_extruded_3d_mesh<I, F>::element_for_interval(int d, int t, std::function<void(I)> f, 
    int xl, int nthreads, bool affinity) const
{
  parallel_for(n_interval(d), [&](int i) {
    const auto id = i + n_ordinal(d) + t * n(d);
    f(id);
  }, xl, nthreads, affinity);

  // for (auto i = 0; i < n_interval(d); i ++) {
  //   const auto id = i + n_ordinal(d) + t * n(d);
  //   f(id);
    // if (d == 2) side_of(2, id);
  // }
}

template <typename I, typename F>
void simplicial_unstructured_extruded_3d_mesh<I, F>::get_simplex(int d, I k, I verts[]) const
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
void simplicial_unstructured_extruded_3d_mesh<I, F>::get_coords(I i, F coords[]) const
{
  const I k = flat_vertex_id(i),
          t = flat_vertex_time(i);
  m.get_coords(k, coords);
  coords[2] = t;
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_extruded_3d_mesh<I, F>::sides(int d, I k) const 
{
  const I i = mod(k, n(d)), t = std::floor((double)k / n(d));
  std::set<I> results;

  I v[5];
  get_simplex(d, i, v);

  if (d == 4) { // currently only pentachrora are supported
    const int type = i / m.n(3);
    if (type == 0) { // 0 1 2 3 3'
      // 0 1 2 3
      // 0 1 2 3'
      { // tet 0123 in original mesh
        const I ot[4] = {v[0], v[1], v[2], v[3]}; // tet in the original mesh
        I otid;
        bool found = m.find_tetrahedra(ot, otid);
        assert(found);

        results.insert(otid + t*n(3)); // 0 1 2 3
        results.insert(otid + t*n(3) + m.n(3)); // 0 1 2 3'
      }
      // 0 1 3 3'
      { // triangle 013 in original mesh
        const I ot[3] = {v[0], v[1], v[3]};
        I otid;
        bool found = m.find_triangle(ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3)); // 0 1 2 2' type
      }
      // 0 2 3 3'
      { // triangle 023 in original mesh
        const I ot[3] = {v[0], v[2], v[3]};
        I otid;
        bool found = m.find_triangle(ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3)); // 0 1 2 2' type
      }
      // 1 2 3 3'
      { // triangle 123 in original mesh
        const I ot[3] = {v[1], v[2], v[3]};
        I otid;
        bool found = m.find_triangle(ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3)); // 0 1 2 2' type
      }
    } else if (type == 1) { // 0 1 2 2'3'
      // 0 1 2 2'
      { // tri 01
        const I ot[3] = {v[0], v[1], v[2]};
        I otid;
        bool found = m.find_triangle(ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3)); // 0 1 2 2' type
      }
      // 0 1 2 3'
      // 0 1 2'3'
      { // tet 0123
        const I ot[4] = {v[0], v[1], v[2], mod(v[4], m.n(0))}; // tet in the original mesh
        I otid;
        bool found = m.find_tetrahedra(ot, otid);
        assert(found);

        results.insert(otid + t*n(3) + m.n(3)); // 0 1 2 3'
        results.insert(otid + t*n(3) + 2*m.n(3)); // 0 1 2'3'
      }
      // 0 2 2'3'
      { // tri 023
        const I ot[3] = {v[0], v[2], mod(v[4], m.n(0))};
        I otid;
        bool found = m.find_triangle(ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + m.n(2)); // 0 1 1'2' type
      }
      // 1 2 2'3'
      { // tri 123
        const I ot[3] = {v[1], v[2], mod(v[4], m.n(0))};
        I otid;
        bool found = m.find_triangle(ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + m.n(2)); // 0 1 1'2' type
      }
    } else if (type == 2) { // 0 1 1'2'3'
      // 0 1 1'2'
      { // tri 012
        const I ot[3] = {v[0], v[1], mod(v[3], m.n(0))};
        I otid;
        bool found = m.find_triangle(ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + m.n(2)); // 0 1 1'2' type
      }
      // 0 1 1'3'
      { // tri 013
        const I ot[3] = {v[0], v[1], mod(v[4], m.n(0))};
        I otid;
        bool found = m.find_triangle(ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + m.n(2)); // 0 1 1'2' type
      }
      // 0 1 2'3'
      // 0 1'2'3'
      { // tet 0123
        const I ot[4] = {v[0], v[1], mod(v[3], m.n(0)), mod(v[4], m.n(0))}; // tet in the original mesh
        I otid;
        bool found = m.find_tetrahedra(ot, otid);
        assert(found);

        results.insert(otid + t*n(3) + 2*m.n(3)); // 0 1 2'3'
        results.insert(otid + t*n(3) + 3*m.n(3)); // 0 1'2'3'
      }
      // 1 1'2'3'
      { // tri 123
        const I ot[3] = {v[1], mod(v[3], m.n(0)), mod(v[4], m.n(0))};
        I otid;
        bool found = m.find_triangle(ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + 2*m.n(2)); // 0 0'1'2' type
      }
    } else if (type == 3) { // 0 0'1'2'3'
      // 0 0'1'2'
      { // tri 012
        const I ot[3] = {v[0], mod(v[2], m.n(0)), mod(v[3], m.n(0))};
        I otid;
        bool found = m.find_triangle(ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + 2*m.n(2)); // 0 0'1'2' type
      }
      // 0 0'1'3'
      { // tri 013
        const I ot[3] = {v[0], mod(v[2], m.n(0)), mod(v[4], m.n(0))};
        I otid;
        bool found = m.find_triangle(ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + 2*m.n(2)); // 0 0'1'2' type
      }
      // 0 0'2'3'
      { // tri 023
        const I ot[3] = {v[0], mod(v[3], m.n(0)), mod(v[4], m.n(0))};
        I otid;
        bool found = m.find_triangle(ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + 2*m.n(2)); // 0 0'1'2' type
      }
      // 0 1'2'3'
      // 0'1'2'3'
      { // tet 0123
        const I ot[4] = {v[0], mod(v[2], m.n(0)), mod(v[3], m.n(0)), mod(v[4], m.n(0))}; // tet in the original mesh
        I otid;
        bool found = m.find_tetrahedra(ot, otid);
        assert(found);

        results.insert(otid + t*n(3) + 2*m.n(3)); // 0 1 2'3'
        results.insert(otid + (t+1)*n(3)); // 0'1'2'3'
      }
    }
  } else if (d == 3) {
    const I v0[4] = {mod(v[0], m.n(0)), mod(v[1], m.n(0)), mod(v[2], m.n(0)), mod(v[3], m.n(0))};
    if (i < 4*m.n(3)) { // "tet" type
      I otid;
      bool found = m.find_tetrahedra(v0, otid);
      assert(found);
     
      I tid;
      if (i < m.n(3)) { // 0 1 2 3
        // 0 1 2  (in 0 1 2 type), find tri 012
        {
          const I tri[3] = {v0[0], v0[1], v0[2]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2));
        }
        // 0 1 3  (in 0 1 2 type), find tri 013
        {
          const I tri[3] = {v0[0], v0[1], v0[3]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2));
        }
        // 0 2 3  (in 0 1 2 type), find tri 023
        {
          const I tri[3] = {v0[0], v0[2], v0[3]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2));
        }
        // 1 2 3  (in 0 1 2 type), find tri 123
        {
          const I tri[3] = {v0[1], v0[2], v0[3]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2));
        }
      } else if (i < 2*m.n(3)) { // 0 1 2 3'
        // 0 1 2  (in 0 1 2  type), find tri 012
        {
          const I tri[3] = {v0[0], v0[1], v0[2]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2));
        }
        // 0 1 3' (in 0 1 2' type), find tri 013
        {
          const I tri[3] = {v0[0], v0[1], v0[3]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2) + m.n(2));
        }
        // 0 2 3' (in 0 1 2' type), find tri 023
        {
          const I tri[3] = {v0[0], v0[2], v0[3]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2) + m.n(2));
        }
        // 1 2 3' (in 0 1 2' type), find tri 123
        {
          const I tri[3] = {v0[1], v0[2], v0[3]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2) + m.n(2));
        }
      } else if (i < 3*m.n(3)) { // 0 1 2'3'
        // 0 1 2' (in 0 1 2' type), find tri 012
        {
          const I tri[3] = {v0[0], v0[1], v0[2]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2) + m.n(2));
        }
        // 0 1 3' (in 0 1 2' type), find tri 013
        {
          const I tri[3] = {v0[0], v0[1], v0[3]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2) + m.n(2));
        }
        // 0 2'3' (in 0 1'2' type), find tri 023
        {
          const I tri[3] = {v0[0], v0[2], v0[3]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2) + 2*m.n(2));
        }
        // 1 2'3' (in 0 1'2' type), find tri 123
        {
          const I tri[3] = {v0[1], v0[2], v0[3]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2) + 2*m.n(2));
        }
      } else if (i < 4*m.n(3)) { // 0 1'2'3'
        // 0 1'2' (in 0 1'2' type), find tri 012
        {
          const I tri[3] = {v0[0], v0[1], v0[2]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2) + 2*m.n(2));
        }
        // 0 1'3' (in 0 1'2' type), find tri 013
        {
          const I tri[3] = {v0[0], v0[1], v0[3]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2) + 2*m.n(2));
        }
        // 0 2'3' (in 0 1'2' type), find tri 023
        {
          const I tri[3] = {v0[0], v0[2], v0[3]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + t*n(2) + 2*m.n(2));
        }
        // 1'2'3' (in 0 1 2  type), find tri 123 in t+1
        {
          const I tri[3] = {v0[1], v0[2], v0[3]};
          assert( m.find_triangle(tri, tid) );
          results.insert(tid + (t+1)*n(2));
        }
      }
    } else { // "tri" type
      I vt[3] = {v0[0], v0[1], v0[3]};
      if (vt[0] == vt[1])
        vt[1] = vt[2];

      I otid;
      bool found = m.find_triangle(vt, otid);
      assert(found);

      I oeid;
      if (i < 4*m.n(3) + m.n(2)) { // 0 1 2 2'
        // 0 1 2 , in 0 1 2  type, find tri 012
        results.insert(otid + t*n(2));
        // 0 1 2', in 0 1 2' type, find tri 012
        results.insert(otid + t*n(2) + m.n(2));
        // 0 2 2', in 0 1 1' type, find edge 02
        {
          const I edge[2] = {vt[0], vt[1]};
          assert( m.find_edge(edge, oeid) );
          results.insert(tid + t*n(2) + 3*m.n(2));
        }
        // 1 2 2', in 0 1 1' type, find edge 12
        {
          const I edge[2] = {vt[1], vt[2]};
          assert( m.find_edge(edge, oeid) );
          results.insert(tid + t*n(2) + 3*m.n(2));
        }
      } else if (i < 4*m.n(3) + 2*m.n(2)) { // 0 1 1'2'
        // 0 1 1', in 0 1 1' type, find edge 01
        {
          const I edge[2] = {vt[0], vt[1]};
          assert( m.find_edge(edge, oeid) );
          results.insert(tid + t*n(2) + 3*m.n(2));
        }
        // 0 1 2', in 0 1 2' type, find tri 012
        results.insert(otid + t*n(2) + m.n(2));
        // 0 1'2', in 0 1'2' type, find tri 012
        results.insert(otid + t*n(2) + 2*m.n(2));
        // 1 1'2', in 0 0'1' type, find edge 12
        {
          const I edge[2] = {vt[1], vt[2]};
          assert( m.find_edge(edge, oeid) );
          results.insert(tid + t*n(2) + 3*m.n(2) + m.n(1));
        }
      } else { // 0 1 1'2'
        // 0 1 1', in 0 1 1' type, find edge 01
        {
          const I edge[2] = {vt[0], vt[1]};
          assert( m.find_edge(edge, oeid) );
          results.insert(tid + t*n(2) + 3*m.n(2));
        }
        // 0 1 2', in 0 1 2' type, find tri 012
        results.insert(otid + t*n(2) + m.n(2));
        // 0 1'2', in 0 1'2' type, find tri 012
        results.insert(otid + t*n(2) + 2*m.n(2));
        // 1 1'2', in 0 0'1' type, find edge 12
        {
          const I edge[2] = {vt[1], vt[2]};
          assert( m.find_edge(edge, oeid) );
          results.insert(tid + t*n(2) + 3*m.n(2) + m.n(1));
        }
      }
    }
  }

  return results;
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_extruded_3d_mesh<I, F>::side_of(int d, I k) const 
{
  const I i = mod(k, n(d)), t = std::floor((double)k / n(d));
  // fprintf(stderr, "k=%d, i=%d, t=%d, n(d)=%d\n", k, i, t, n(d));
  std::set<I> results;

  I v[5];
  get_simplex(d, i, v);

  if (d == 3) {
    const I v0[4] = {mod(v[0], m.n(0)), mod(v[1], m.n(0)), mod(v[2], m.n(0)), mod(v[3], m.n(0))};
    if (i < 4*m.n(3)) { // "tet" type
      I otid;
      bool found = m.find_tetrahedra(v0, otid);
      assert(found);

      if (i < m.n(3)) { // 0 1 2 3
        // t:   0 1 2 3 3' 
        results.insert(otid + t*n(4));
        // t-1: 0 0'1'2'3'
        results.insert(otid + (t-1)*n(4) + 3*m.n(3));
      } else if (i < 2*m.n(3)) { // 0 1 2 3'
        // t:   0 1 2 3 3'
        results.insert(otid + t*n(4));
        // t:   0 1 2 2'3'
        results.insert(otid + t*n(4) + m.n(3));
      } else if (i < 3*m.n(3)) { // 0 1 2'3'
        // t:   0 1 2 2'3'
        results.insert(otid + t*n(4) + m.n(3));
        // t:   0 1 1'2'3'
        results.insert(otid + t*n(4) + 2*m.n(3));
      } else if (i < 4*m.n(3)) { // 0 1'2'3'
        // t:   0 1 1'2'3'
        results.insert(otid + t*n(4) + 2*m.n(3));
        // t:   0 0'1'2'3'
        results.insert(otid + t*n(4) + 3*m.n(3));
      }
    } else { // "tri" type
      I vt[3] = {mod(v[0], m.n(0)), mod(v[1], m.n(0)), mod(v[3], m.n(0))};
      if (vt[0] == vt[1])
        vt[1] = mod(v[2], m.n(0));

      I otid;
      bool found = m.find_triangle(vt, otid);
      assert(found);

      // TODO: also consider "pos"
      if (i < 4*m.n(3) + m.n(2)) { // 0 1 2 2'
        // find which tet in the base mesh has the 012 face
        // 0 1 2 2'3'
        // 0 1 2 3 3' in another prism
      } else if (i < 4*m.n(3) + 2*m.n(2)) { // 0 1 1'2'
        // 0 1 1'2'3'
        // 0 1 2 2'3' in another prism
      } else { // 0 0'1'2'
        // 0 0'1'2'3'
        // 0 1 1'2'3' in another prism
      }
    }
  } else if (d == 2) {
    // find all two tets in the base mesh that contains tri 012
    // check "position" too, e,g, 013 in 0123
    { // type 0 1 2 (or 0'1'2'), pos=0, 1, 2, 3
      // 0 1 2 3
      // 0 1 2 3'
      // 0 1 2 2'
      // 0 0'1'2', t-1
    }
    { // type 0 1 2', pos=0, 1, 2, 3
      // 0 1 2'3'
      // 0 1 2 2'
    }
    { // type 0 1'2', pos=0, 1, 2, 3
      // 0 1 1'2'
      // 0 1'2'3'
    }
    // find all triangles in the base mesh that contains edge 01
    { // type 0 1 1', pos=0, 1, 2
      // 0 1 1'2'
      // 0 0'1'2'
    }
    { // type 0 0'1', pos=0, 1, 2
      // 0 0'1'2'
    }
  }

  return results;
}

template <typename I, typename F>
int simplicial_unstructured_extruded_3d_mesh<I, F>::face_type(I k) const
{
  const I i = mod(k, n(2));

  if (i < m.n(2)) return 0;  // I
  else if (i < 2*m.n(2)) return 1; // II
  else if (i < 3*m.n(2)) return 2; // III 
  else if (i < 3*m.n(2) + m.n(1)) return 3; // IV
  else if (i < 3*m.n(2) + 2*m.n(1)) return 4; // V
  else return -1;
}
#endif

}

#endif
