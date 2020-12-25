#ifndef _HYPERMESH_simplicial_unstructured_extruded_3d_mesh_HH
#define _HYPERMESH_simplicial_unstructured_extruded_3d_mesh_HH

#include <ftk/mesh/simplicial_unstructured_3d_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_unstructured_extruded_3d_mesh : public object { // extruded from 
  simplicial_unstructured_extruded_3d_mesh(const simplicial_unstructured_3d_mesh<I, F>& m_) : m(m_) {extrude();}

  size_t n(int d) const;
  size_t n_ordinal(int d) const { if (d == 4) return m.n(3)*4; else return m.n(d); }
  size_t n_interval(int d) const;

  bool is_ordinal(int d, I i) const {
    I k = mod(i, n(d));
    return k < n_ordinal(d);
  }

  I flat_vertex_id(I i) const { return mod(i, m.n(0)); }
  I flat_vertex_time(I i) const { return std::floor((double)i / n(0)); } // return i / m.n(0); }
  I extruded_vertex_id(I i, bool t=true) { return t ? i + m.n(0) : i; }

  int tet_type(I i) const { return simplex_type(3, i); }
  int simplex_type(int d, I i) const;

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

#if 0
  ftk::ndarray<I> pents,// {5, n(4)}
                  tets, // {4, n(3)}
                  tris, // {3, n(2)}
                  edges;// {2, n(1)}

  ftk::ndarray<I> pents_sides,// {5, n(4)}
                  tets_sides, // {4, n(3)}
                  tris_sides; // {3, n(2)}
#endif
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
  fprintf(stderr, "4d mesh initialized, #pent=%zu, #tet=%zu, #tri=%zu, #edge=%zu, #vert=%zu\n", 
      n(4), n(3), n(2), n(1), n(0));
#if 0
  // pents
  pents.reshape({5, n(4)});
  for (auto i = 0; i < m.n(3); i ++) {
    I tet[4];
    m.get_simplex(3, i, tet);

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
    m.get_simplex(3, i, tet);

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
    m.get_simplex(2, i, tri);

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
    m.get_simplex(2, i, tri);

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
    m.get_simplex(1, i, edge);

    tris(0, i+3*m.n(2)) = edge[0]; // 0 1 1'
    tris(1, i+3*m.n(2)) = edge[1];
    tris(2, i+3*m.n(2)) = edge[1] + m.n(0);

    tris(0, i+3*m.n(2)+m.n(1)) = edge[0]; // 0 0'1' // tri_type==4
    tris(1, i+3*m.n(2)+m.n(1)) = edge[0] + m.n(0);
    tris(2, i+3*m.n(2)+m.n(1)) = edge[1] + m.n(0);
  }

  // unique edges
  edges.reshape({2, n(1)});
  for (auto i = 0; i < m.n(1); i ++) {
    I v[2];
    m.get_simplex(1, i, v);

    edges(0, i) = v[0]; // base edge
    edges(1, i) = v[1];

    edges(0, i+m.n(1)) = v[0]; // 0 1'
    edges(1, i+m.n(1)) = v[1] + m.n(0);
  }
  for (auto i = 0; i < m.n(0); i ++) {
    edges(0, i+2*m.n(1)) = i;
    edges(1, i+2*m.n(1)) = i + m.n(0);
  }
#endif
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
  const int type = simplex_type(d, i);
  
  // fprintf(stderr, "d=%d, k=%d, i=%d, t=%d, type=%d, offset=%d\n", d, k, i, t, type, offset);
  if (d == 0) {
    verts[0] = i; // + offset;
  } else if (d == 1) {
    if (type < 2) {
      I edge[2];
      m.get_simplex(1, i % m.n(1), edge);

      switch (type) {
        case 0: // 0 1
          verts[0] = edge[0];
          verts[1] = edge[1];
          break;

        case 1: // 0 1'
          verts[0] = edge[0];
          verts[1] = edge[1] + m.n(0);
      }
    } else {
      verts[0] = i - 2*m.n(1);
      verts[1] = i - 2*m.n(1) + m.n(0);
    }
  } else if (d == 2) {
    if (type < 3) {
      I tri[3];
      m.get_simplex(2, i % m.n(2), tri);
      switch (type) {
        case 0: // 0 1 2
          verts[0] = tri[0];
          verts[1] = tri[1];
          verts[2] = tri[2];
          break;

        case 1: // 0 1 2'
          verts[0] = tri[0];
          verts[1] = tri[1];
          verts[2] = tri[2] + m.n(0);
          break;

        case 2: // 0 1'2'
          verts[0] = tri[0];
          verts[1] = tri[1] + m.n(0);
          verts[2] = tri[2] + m.n(0);
          break;

        default: assert(false);
      }
    } else {
      I edge[2];
      m.get_simplex(1, (i - 3*m.n(2)) % m.n(1), edge);
      switch (type) {
        case 3: 
          verts[0] = edge[0];
          verts[1] = edge[1];
          verts[2] = edge[1] + m.n(0);
          break;

        case 4:
          verts[0] = edge[0];
          verts[1] = edge[0] + m.n(0);
          verts[2] = edge[1] + m.n(0);
          break;

        default: assert(false);
      }
    }
  } else if (d == 3) {
    if (type < 4) {
      I tet[4];
      m.get_simplex(3, i % m.n(3), tet);
      switch (type) {
        case 0: // 0 1 2 3, "base"
          verts[0] = tet[0];
          verts[1] = tet[1];
          verts[2] = tet[2];
          verts[3] = tet[3];
          break;

        case 1: // 0 1 2 3'
          verts[0] = tet[0];
          verts[1] = tet[1];
          verts[2] = tet[2];
          verts[3] = tet[3] + m.n(0);
          break;

        case 2: // 0 1 2'3'
          verts[0] = tet[0];
          verts[1] = tet[1];
          verts[2] = tet[2] + m.n(0);
          verts[3] = tet[3] + m.n(0);
          break;

        case 3: // 0 1'2'3'
          verts[0] = tet[0];
          verts[1] = tet[1] + m.n(0);
          verts[2] = tet[2] + m.n(0);
          verts[3] = tet[3] + m.n(0);
          break;

        default: 
          assert(false);
      }
    } else {
      I tri[3];
      m.get_simplex(2, (i - 4*m.n(3)) % m.n(2), tri);
      switch (type) {
        case 4: // 0 1 2 2'
          verts[0] = tri[0];
          verts[1] = tri[1];
          verts[2] = tri[2];
          verts[3] = tri[2] + m.n(0);
          break;

        case 5: // 0 1 1'2'
          verts[0] = tri[0];
          verts[1] = tri[1];
          verts[2] = tri[1] + m.n(0);
          verts[3] = tri[2] + m.n(0);
          break;

        case 6: // 0 0'1'2'
          verts[0] = tri[0];
          verts[1] = tri[0] + m.n(0);
          verts[2] = tri[1] + m.n(0);
          verts[3] = tri[2] + m.n(0);
          break;

        default:
          assert(false);
      }
    }
  } else if (d == 4) {
    I tet[4];
    m.get_simplex(3, i % m.n(3), tet);

    switch (type) {
      case 0: // type I:   0 1 2 3 3'
        verts[0] = tet[0]; 
        verts[1] = tet[1];
        verts[2] = tet[2];
        verts[3] = tet[3];
        verts[4] = tet[3] + m.n(0);
        break;

      case 1: // type II:  0 1 2 2'3'
        verts[0] = tet[0];
        verts[1] = tet[1];
        verts[2] = tet[2];
        verts[3] = tet[2] + m.n(0);
        verts[4] = tet[3] + m.n(0);
        break;

      case 2: // type III: 0 1 1'2'3'
        verts[0] = tet[0];
        verts[1] = tet[1];
        verts[2] = tet[1] + m.n(0);
        verts[3] = tet[2] + m.n(0);
        verts[4] = tet[3] + m.n(0);
        break;

      case 3: // type IV:  0 1'1'2'3'
        verts[0] = tet[0];
        verts[1] = tet[1] + m.n(0);
        verts[2] = tet[1] + m.n(0);
        verts[3] = tet[2] + m.n(0);
        verts[4] = tet[3] + m.n(0);
        break;

      default: 
        assert(false);
    }
  }

  for (int i = 0; i < d+1; i ++)
    verts[i] += offset;
}

template <typename I, typename F>
void simplicial_unstructured_extruded_3d_mesh<I, F>::get_coords(I i, F coords[]) const
{
  const I k = flat_vertex_id(i),
          t = flat_vertex_time(i);
  m.get_coords(k, coords);
  coords[3] = t;
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_extruded_3d_mesh<I, F>::sides(int d, I k) const 
{
  const I i = mod(k, n(d)), t = std::floor((double)k / n(d));
  std::set<I> results;

  I v[5];
  get_simplex(d, i, v);

  if (d == 4) { // currently only pentachrora and tetrahedra are supported
    const int type = i / m.n(3);
    if (type == 0) { // 0 1 2 3 3'
      // 0 1 2 3
      // 0 1 2 3'
      { // tet 0123 in original mesh
        const I ot[4] = {v[0], v[1], v[2], v[3]}; // tet in the original mesh
        I otid;
        bool found = m.find_simplex(3, ot, otid);
        assert(found);

        results.insert(otid + t*n(3)); // 0 1 2 3
        results.insert(otid + t*n(3) + m.n(3)); // 0 1 2 3'
      }
      // 0 1 3 3'
      { // triangle 013 in original mesh
        const I ot[3] = {v[0], v[1], v[3]};
        I otid;
        bool found = m.find_simplex(2, ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3)); // 0 1 2 2' type
      }
      // 0 2 3 3'
      { // triangle 023 in original mesh
        const I ot[3] = {v[0], v[2], v[3]};
        I otid;
        bool found = m.find_simplex(2, ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3)); // 0 1 2 2' type
      }
      // 1 2 3 3'
      { // triangle 123 in original mesh
        const I ot[3] = {v[1], v[2], v[3]};
        I otid;
        bool found = m.find_simplex(2, ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3)); // 0 1 2 2' type
      }
    } else if (type == 1) { // 0 1 2 2'3'
      // 0 1 2 2'
      { // tri 012
        const I ot[3] = {v[0], v[1], v[2]};
        I otid;
        bool found = m.find_simplex(2, ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3)); // 0 1 2 2' type
      }
      // 0 1 2 3'
      // 0 1 2'3'
      { // tet 0123
        const I ot[4] = {v[0], v[1], v[2], mod(v[4], m.n(0))}; // tet in the original mesh
        I otid;
        bool found = m.find_simplex(3, ot, otid);
        assert(found);

        results.insert(otid + t*n(3) + m.n(3)); // 0 1 2 3'
        results.insert(otid + t*n(3) + 2*m.n(3)); // 0 1 2'3'
      }
      // 0 2 2'3'
      { // tri 023
        const I ot[3] = {v[0], v[2], mod(v[4], m.n(0))};
        I otid;
        bool found = m.find_simplex(2, ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + m.n(2)); // 0 1 1'2' type
      }
      // 1 2 2'3'
      { // tri 123
        const I ot[3] = {v[1], v[2], mod(v[4], m.n(0))};
        I otid;
        bool found = m.find_simplex(2, ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + m.n(2)); // 0 1 1'2' type
      }
    } else if (type == 2) { // 0 1 1'2'3'
      // 0 1 1'2'
      { // tri 012
        const I ot[3] = {v[0], v[1], mod(v[3], m.n(0))};
        I otid;
        bool found = m.find_simplex(2, ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + m.n(2)); // 0 1 1'2' type
      }
      // 0 1 1'3'
      { // tri 013 
        const I ot[3] = {v[0], v[1], mod(v[4], m.n(0))};
        I otid;
        bool found = m.find_simplex(2, ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + m.n(2)); // 0 1 1'2' type
      }
      // 0 1 2'3'
      // 0 1'2'3'
      { // tet 0123
        const I ot[4] = {v[0], v[1], mod(v[3], m.n(0)), mod(v[4], m.n(0))}; // tet in the original mesh
        I otid;
        bool found = m.find_simplex(3, ot, otid);
        assert(found);

        results.insert(otid + t*n(3) + 2*m.n(3)); // 0 1 2'3'
        results.insert(otid + t*n(3) + 3*m.n(3)); // 0 1'2'3'
      }
      // 1 1'2'3'
      { // tri 123 
        const I ot[3] = {v[1], mod(v[3], m.n(0)), mod(v[4], m.n(0))};
        I otid;
        bool found = m.find_simplex(2, ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + 2*m.n(2)); // 0 0'1'2' type
      }
    } else if (type == 3) { // 0 0'1'2'3'
      // 0 0'1'2'
      { // tri 012
        const I ot[3] = {v[0], mod(v[2], m.n(0)), mod(v[3], m.n(0))};
        I otid;
        bool found = m.find_simplex(2, ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + 2*m.n(2)); // 0 0'1'2' type
      }
      // 0 0'1'3'
      { // tri 013
        const I ot[3] = {v[0], mod(v[2], m.n(0)), mod(v[4], m.n(0))};
        I otid;
        bool found = m.find_simplex(2, ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + 2*m.n(2)); // 0 0'1'2' type
      }
      // 0 0'2'3'
      { // tri 023
        const I ot[3] = {v[0], mod(v[3], m.n(0)), mod(v[4], m.n(0))};
        I otid;
        bool found = m.find_simplex(2, ot, otid);
        assert(found);
        results.insert(otid + t*n(3) + 4*m.n(3) + 2*m.n(2)); // 0 0'1'2' type
      }
      // 0 1'2'3'
      // 0'1'2'3'
      { // tet 0123
        const I ot[4] = {v[0], mod(v[2], m.n(0)), mod(v[3], m.n(0)), mod(v[4], m.n(0))}; // tet in the original mesh
        I otid;
        bool found = m.find_simplex(3, ot, otid);
        assert(found);

        results.insert(otid + t*n(3) + 3*m.n(3)); // 0 1'2'3'
        results.insert(otid + (t+1)*n(3)); // 0'1'2'3'
      }
    }
  } else if (d == 3) {
    I v0[4] = {mod(v[0], m.n(0)), mod(v[1], m.n(0)), mod(v[2], m.n(0)), mod(v[3], m.n(0))};
    if (i < 4*m.n(3)) { // "tet" type
      I otid;
      bool found = m.find_simplex(3, v0, otid);
      assert(found);
     
      I tid;
      if (i < m.n(3)) { // 0 1 2 3
        // 0 1 2  (in 0 1 2 type), find tri 012
        {
          const I tri[3] = {v0[0], v0[1], v0[2]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2));
        }
        // 0 1 3  (in 0 1 2 type), find tri 013
        {
          const I tri[3] = {v0[0], v0[1], v0[3]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2));
        }
        // 0 2 3  (in 0 1 2 type), find tri 023
        {
          const I tri[3] = {v0[0], v0[2], v0[3]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2));
        }
        // 1 2 3  (in 0 1 2 type), find tri 123
        {
          const I tri[3] = {v0[1], v0[2], v0[3]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2));
        }
      } else if (i < 2*m.n(3)) { // 0 1 2 3'
        // 0 1 2  (in 0 1 2  type), find tri 012
        {
          const I tri[3] = {v0[0], v0[1], v0[2]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2));
        }
        // 0 1 3' (in 0 1 2' type), find tri 013
        {
          const I tri[3] = {v0[0], v0[1], v0[3]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2) + m.n(2));
        }
        // 0 2 3' (in 0 1 2' type), find tri 023
        {
          const I tri[3] = {v0[0], v0[2], v0[3]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2) + m.n(2));
        }
        // 1 2 3' (in 0 1 2' type), find tri 123
        {
          const I tri[3] = {v0[1], v0[2], v0[3]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2) + m.n(2));
        }
      } else if (i < 3*m.n(3)) { // 0 1 2'3' // tet_type==2
        // 0 1 2' (in 0 1 2' type), find tri 012
        {
          const I tri[3] = {v0[0], v0[1], v0[2]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2) + m.n(2));
        }
        // 0 1 3' (in 0 1 2' type), find tri 013
        {
          const I tri[3] = {v0[0], v0[1], v0[3]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2) + m.n(2));
        }
        // 0 2'3' (in 0 1'2' type), find tri 023  
        {
          const I tri[3] = {v0[0], v0[2], v0[3]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2) + 2*m.n(2));
        }
        // 1 2'3' (in 0 1'2' type), find tri 123  
        {
          const I tri[3] = {v0[1], v0[2], v0[3]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2) + 2*m.n(2));
        }
      } else if (i < 4*m.n(3)) { // 0 1'2'3' // tet_type==3
        // 0 1'2' (in 0 1'2' type), find tri 012
        {
          const I tri[3] = {v0[0], v0[1], v0[2]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2) + 2*m.n(2));
        }
        // 0 1'3' (in 0 1'2' type), find tri 013
        {
          const I tri[3] = {v0[0], v0[1], v0[3]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2) + 2*m.n(2));
        }
        // 0 2'3' (in 0 1'2' type), find tri 023
        {
          const I tri[3] = {v0[0], v0[2], v0[3]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + t*n(2) + 2*m.n(2));
        }
        // 1'2'3' (in 0 1 2  type), find tri 123 in t+1
        {
          const I tri[3] = {v0[1], v0[2], v0[3]};
          bool found = m.find_simplex(2, tri, tid);
          assert(found);
          results.insert(tid + (t+1)*n(2));
        }
      }
    } else { // "tri" type
      // I vt[3] = {v0[0], v0[1], v0[3]};
      // if (vt[0] == vt[1])
      //   vt[1] = vt[2];

      I otid, oeid;
      if (i < 4*m.n(3) + m.n(2)) { // 0 1 2 2'
        I vt[3] = {v0[0], v0[1], v0[2]};
        bool found = m.find_simplex(2, vt, otid);
        assert(found);

        // 0 1 2 , in 0 1 2  type, find tri 012
        results.insert(otid + t*n(2));
        // 0 1 2', in 0 1 2' type, find tri 012
        results.insert(otid + t*n(2) + m.n(2));
        // 0 2 2', in 0 1 1' type, find edge 02
        {
          const I edge[2] = {v0[0], v0[2]};
          bool found = m.find_simplex(1, edge, oeid);
          assert(found);
          results.insert(oeid + t*n(2) + 3*m.n(2));
        }
        // 1 2 2', in 0 1 1' type, find edge 12
        {
          const I edge[2] = {v0[1], v0[2]};
          bool found = m.find_simplex(1, edge, oeid);
          assert(found);
          results.insert(oeid + t*n(2) + 3*m.n(2));
        }
      } else if (i < 4*m.n(3) + 2*m.n(2)) { // 0 1 1'2' // tet_type==5
        I vt[3] = {v0[0], v0[1], v0[3]};
        bool found = m.find_simplex(2, vt, otid);
        assert(found);
        
        // 0 1 1', in 0 1 1' type, find edge 01
        {
          const I edge[2] = {v0[0], v0[1]};
          bool found = m.find_simplex(1, edge, oeid);
          assert(found);
          results.insert(oeid + t*n(2) + 3*m.n(2));
        }
        // 0 1 2', in 0 1 2' type, find tri 012
        results.insert(otid + t*n(2) + m.n(2));
        // 0 1'2', in 0 1'2' type, find tri 012
        results.insert(otid + t*n(2) + 2*m.n(2)); 
        // 1 1'2', in 0 0'1' type, find edge 12
        {
          const I edge[2] = {v0[1], v0[3]};
          bool found = m.find_simplex(1, edge, oeid);
          assert(found);
          results.insert(oeid + t*n(2) + 3*m.n(2) + m.n(1));
        }
      } else { // 0 0'1'2' // tet_type==6
        I vt[3] = {v0[0], v0[2], v0[3]};
        bool found = m.find_simplex(2, vt, otid);
        if (!found) {
          fprintf(stderr, "FATAL:, v0=%d, %d, %d, %d, vt=%d, %d, %d\n", 
              v0[0], v0[1], v0[2], v0[3], 
              vt[0], vt[1], vt[2]);
        }
        assert(found);
      
        // 0 0'1', in 0 0'1' type, find edge 01
        {
          const I edge[2] = {v0[0], v0[2]};
          bool found = m.find_simplex(1, edge, oeid);
          assert(found);
          results.insert(oeid + t*n(2) + 3*m.n(2) + m.n(1));
        }
        // 0 0'2', in 0 0'1' type, find edge 02
        {
          const I edge[2] = {v0[0], v0[3]};
          bool found = m.find_simplex(1, edge, oeid);
          assert(found);
          results.insert(oeid + t*n(2) + 3*m.n(2) + m.n(1));
        }
        // 0 1'2', in 0 1'2' type, find tri 012
        results.insert(otid + t*n(2) + 2*m.n(2));
        // 0'1'2', in 0 1 2  type in t+1, find tr 012
        results.insert(otid + (t+1)*n(2));
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
  
  auto position12 = [](const I edge[2], const I tri[3]) { // edgie is a side of tri
    for (int i = 0; i < 3; i ++) {
      bool f = false;
      for (int j = 0; j < 2; j ++)
        if (tri[i] == edge[j]) {
          f = true;
          break;
        }
      if (!f) return i;
    }
    assert(false);
    return -1;
  };
      
  auto position23 = [](const I tri[3], const I tet[4]) { // tri is a side of tet
    // fprintf(stderr, "TRI=%d, %d, %d, TET=%d, %d, %d, %d\n", tri[0], tri[1], tri[2], tet[0], tet[1], tet[2], tet[3]);
    for (int i = 0; i < 4; i ++) {
      bool f = false;
      for (int j = 0; j < 3; j ++)
        if (tet[i] == tri[j]) {
          f = true;
          break;
        }
      if (!f) return i;
    }
    assert(false);
    return -1;
  };
      

  if (d == 3) {
    const I v0[4] = {mod(v[0], m.n(0)), mod(v[1], m.n(0)), mod(v[2], m.n(0)), mod(v[3], m.n(0))};
    // fprintf(stderr, "INPUT TET: v=%d, %d, %d, %d, v0=%d, %d, %d, %d, type=%d\n", 
    //     v[0], v[1], v[2], v[3],
    //     v0[0], v0[1], v0[2], v0[3], tet_type(i));
    if (i < 4*m.n(3)) { // "tet" type
      I otid;
      bool found = m.find_simplex(3, v0, otid);
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

        // fprintf(stderr, "%d, %d\n", 
        //     otid + t*n(4) + m.n(3),
        //     otid + t*n(4) + 2*m.n(3));

      } else if (i < 4*m.n(3)) { // 0 1'2'3'
        // t:   0 1 1'2'3'
        results.insert(otid + t*n(4) + 2*m.n(3));
        // t:   0 0'1'2'3'
        results.insert(otid + t*n(4) + 3*m.n(3));
      }
    } else { // "tri" type
      // fprintf(stderr, "TRI TYPE!!!\n");
      I otid; // triangle id in the base mesh
     
      if (i < 4*m.n(3) + m.n(2)) { // 0 1 2 2'
        const I vt[3] = {v0[0], v0[1], v0[2]};
        bool found = m.find_simplex(2, vt, otid);
        assert(found);

        for (const auto &otetid : m.side_of(2, otid)) {
          I tet[4];
          m.get_simplex(3, otetid, tet);
          int pos = position23(vt, tet);

          // fprintf(stderr, "POS=%d\n", pos);

          if (pos == 0) 
            results.insert(otetid + t*n(4));
          else if (pos == 1)
            results.insert(otetid + t*n(4));
          else if (pos == 2)
            results.insert(otetid + t*n(4));
          else if (pos == 3)
            results.insert(otetid + t*n(4) + m.n(3));

#if 0
          bool upper = (vt[0] == tet[0]);
          fprintf(stderr, "OTETID=%d, pos=%d, upper=%d, tet=%d, %d, %d, %d\n", 
              otetid, pos, upper, tet[0], tet[1], tet[2], tet[3]);
          if (upper) {
            results.insert(otetid + t*n(4) + m.n(3));
            // fprintf(stderr, "adding upper: %d\n", otetid + t*n(4)); //  + m.n(3));
          }
          else { // lower
            results.insert(otetid + t*n(4)); //  + m.n(3));
            // fprintf(stderr, "adding lower: %d\n", otetid + t*n(4));
          }
#endif
          // results.insert(otetid + t*n(4));
        }
        // find which tet in the base mesh has the 012 face
        // 0 1 2 2'3'
        // 0 1 2 3 3' in another prism
      } else if (i < 4*m.n(3) + 2*m.n(2)) { // 0 1 1'2'
        const I vt[3] = {v0[0], v0[1], v0[3]};
        bool found = m.find_simplex(2, vt, otid);
        assert(found);
        // 0 1 1'2'3'
        // 0 1 2 2'3' in another prism
        for (const auto &otetid : m.side_of(2, otid)) {
          I tet[4];
          m.get_simplex(3, otetid, tet);
          
          int pos = position23(vt, tet);
          // fprintf(stderr, "POS=%d\n", pos);

          if (pos == 0) 
            results.insert(otetid + t*n(4) + m.n(3));
          else if (pos == 1)
            results.insert(otetid + t*n(4) + m.n(3));
          else if (pos == 2)
            results.insert(otetid + t*n(4) + 2*m.n(3));
          else if (pos == 3)
            results.insert(otetid + t*n(4) + 2*m.n(3));
        }
      } else { // 0 0'1'2'
        const I vt[3] = {v0[0], v0[2], v0[3]};
        bool found = m.find_simplex(2, vt, otid);
        assert(found);
        // 0 0'1'2'3'
        // 0 1 1'2'3' in another prism
        for (const auto &otetid : m.side_of(2, otid)) {
          I tet[4];
          m.get_simplex(3, otetid, tet);
          
          int pos = position23(vt, tet);
          // fprintf(stderr, "POS=%d\n", pos);

          if (pos == 0) 
            results.insert(otetid + t*n(4) + 2*m.n(3));
          else if (pos == 1)
            results.insert(otetid + t*n(4) + 3*m.n(3));
          else if (pos == 2)
            results.insert(otetid + t*n(4) + 3*m.n(3));
          else if (pos == 3)
            results.insert(otetid + t*n(4) + 3*m.n(3));
#if 0
          bool upper = (vt[0] == tet[0]);
          if (upper)
            results.insert(otetid + t*n(4) + 3*m.n(3));
          else // lower
            results.insert(otetid + t*n(4) + 2*m.n(3));
#endif
        }
      }
    }
  } else if (d == 2) { // WIP
    // find all tets in the base mesh that contains tri 012
    // check "position" too, e,g, 013 in 0123
    const I v0[3] = {mod(v[0], m.n(0)), mod(v[1], m.n(0)), mod(v[2], m.n(0))};
    I otid, oeid;

    if (i < m.n(2)) { // type 0 1 2 (or 0'1'2'), pos=0, 1, 2, 3
      bool found = m.find_simplex(2, v0, otid);
      assert(found);

      for (auto otetid : m.side_of(2, otid)) {
        I tet[4];
        m.get_simplex(3, otetid, tet);

        int pos = position23(v0, tet);
        if (pos == 0) { // a 0 1 2 in base mesh, find all tets that contains 012 or 0'1'2' in the a012--a'0'1'2' prism
                        // equiv to find 123 in 0123-0'1'2'3'
          // 0 1 2 3
          results.insert(otetid + t*n(3));
          // 0 1'2'3' (t-1)
          results.insert(otetid + (t-1)*n(3) + 3*m.n(3));
        } else if (pos == 1) { // 0 a 1 2, equil. 023 in 0123-0'1'2'3'
          // 0 1 2 3
          results.insert(otetid + t*n(3));
          // 0 2 3 3' (0 1 2 2' type)
          {
            int otriid;
            I otri[3] = {tet[0], tet[2], tet[3]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3));
          }
        } else if (pos == 2) { // 0 1 a 2, equil. 013 in 0123-0'1'2'3'
          // 0 1 2 3
          results.insert(otetid + t*n(3));
          // 0 1 3 3' (0 1 2 2' type)
          {
            int otriid;
            I otri[3] = {tet[0], tet[1], tet[3]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3));
          }
        } else if (pos == 3) { // 0 1 2 a, equil. 012 in 0123-0'1'2'3'
          // 0 1 2 3
          results.insert(otetid + t*n(3));
          // 0 1 2 2' (0 1 2 2' type)
          {
            int otriid;
            I otri[3] = {tet[0], tet[1], tet[2]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3));
          }
          // 0 1 2 3'
          results.insert(otetid + t*n(3) + m.n(3));
          // 0 0'1'2' (t-1, 0 0'1'2' type)
          {
            int otriid;
            I otri[3] = {tet[0], tet[1], tet[2]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + (t-1)*n(3) + 4*m.n(3) + 2*m.n(2));
          }
        }
      }
    } else if (i < 2*m.n(2)) { // type 0 1 2', pos=0, 1, 2, 3
      bool found = m.find_simplex(2, v0, otid);
      assert(found);

      for (auto otetid : m.side_of(2, otid)) {
        I tet[4];
        m.get_simplex(3, otetid, tet);

        int pos = position23(v0, tet);
        if (pos == 0) { // a 0 1 2 in base mesh, find all tets that contains 012' in the a012--a'0'1'2' prism
                        // equiv 123' in 0123-0'1'2'3'
          // 0 1 2 3'
          results.insert(otetid + t*n(3) + m.n(3));
          // 1 2 2'3' (0 1 1'2' type)
          {
            int otriid;
            I otri[3] = {tet[1], tet[2], tet[3]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3) + m.n(2));
          }
        } else if (pos == 1) { // 0 a 1 2, find 012', equiv to finding 023' in 0123-0'1'2'3'
          // 0 1 2 3'
          results.insert(otetid + t*n(3) + m.n(3));
          // 0 2 2'3' (0 1 1'2' type)
          {
            int otriid;
            I otri[3] = {tet[0], tet[2], tet[3]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3) + m.n(2));
          }
        } else if (pos == 2) { // 0 1 a 2, equiv to finding 013' in 0123-0'1'2'3'
          // 0 1 1'3' (0 1 1'2' type)
          {
            int otriid;
            I otri[3] = {tet[0], tet[1], tet[3]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3) + m.n(2));
          }
          // 0 1 3 3' (0 1 2 2' type)
          {
            int otriid;
            I otri[3] = {tet[0], tet[1], tet[3]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3));
          }
        } else if (pos == 3) { // 0 1 2 a, equiv to finding 012' in 0123-0'1'2'3'
          // 0 1 2 2' (0 1 2 2' type)
          {
            int otriid;
            I otri[3] = {tet[0], tet[1], tet[2]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3));
          }
          // 0 1 1'2' (0 1 1'2' type)
          {
            int otriid;
            I otri[3] = {tet[0], tet[1], tet[2]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3) + m.n(2));
          }
        }
      }
    } else if (i < 3*m.n(2)) { // type 0 1'2', pos=0, 1, 2, 3
      bool found = m.find_simplex(2, v0, otid);
      assert(found);
      
      for (auto otetid : m.side_of(2, otid)) {
        I tet[4];
        m.get_simplex(3, otetid, tet);

        int pos = position23(v0, tet);
        if (pos == 0) { // 1 2'3' in tet-prism
          // 0 1 2'3'
          results.insert(otetid + t*n(3) + 2*m.n(3));
          // 1 1'2'3' (0 0'1'2' type)
          {
            int otriid;
            I otri[3] = {tet[1], tet[2], tet[3]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3) + 2*m.n(2));
          }
          // 1 2 2'3' (0 1 1'2' type)
          {
            int otriid;
            I otri[3] = {tet[1], tet[2], tet[3]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3) + m.n(2));
          }
        } else if (pos == 1) { // 0 2'3' in tet-prism
          // 0 1 2'3'
          results.insert(otetid + t*n(3) + 2*m.n(3));
          // 0 1'2'3'
          results.insert(otetid + t*n(3) + 3*m.n(3));
          // 0 2 2'3' (0 1 1'2' type)
          {
            int otriid;
            I otri[3] = {tet[0], tet[2], tet[3]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3) + m.n(2));
          }
        } else if (pos == 2) { // 0 1'3' in tet-prism
          // 0 1 1'3' (0 1 1'2' type)
          {
            int otriid;
            I otri[3] = {tet[0], tet[1], tet[3]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3) + m.n(2));
          }
          // 0 0'1'3' (0 0'1'2' type)
          {
            int otriid;
            I otri[3] = {tet[0], tet[1], tet[3]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3) + 2*m.n(2));
          }
          // 0 1'2'3'
          results.insert(otetid + t*n(3) + 3*m.n(3));
        } else if (pos == 3) { // 0 1'2' in tet-prism
          // 0 0'1'2' (0 0'1'2' type)
          {
            int otriid;
            I otri[3] = {tet[0], tet[1], tet[2]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3) + 2*m.n(2));
          }
          // 0 1 1'2' (0 1 1'2' type)
          {
            int otriid;
            I otri[3] = {tet[0], tet[1], tet[2]};
            bool succ = m.find_simplex(2, otri, otriid);
            assert(succ);
            results.insert(otriid + t*n(3) + 4*m.n(3) + m.n(2));
          }
          // 0 1'2'3'
          results.insert(otetid + t*n(3) + 3*m.n(3));
        }
      }
    } else if (i < 3*m.n(2) + m.n(1)) { // find all triangles in the base mesh that contains edge 01 // type 0 1 1', pos=0, 1, 2
      const I ve[2] = {v0[0], v0[1]};
      bool found = m.find_simplex(1, ve, oeid);
      assert(found);

      for (const auto otriid : m.side_of(1, oeid)) {
        I tri[3];
        m.get_simplex(2, otriid, tri);

        const int pos = position12(ve, tri);
        if (pos == 0) { // 1 2 2'
          // 0 1 2 2'
          results.insert(otriid + t*n(3) + 4*m.n(3));
        } else if (pos == 1) { // 0 2 2'
          // 0 1 2 2'
          results.insert(otriid + t*n(3) + 4*m.n(3));
        } else if (pos == 2) { // 0 1 1'
          // 0 1 1'2'
          results.insert(otriid + t*n(3) + 4*m.n(3) + m.n(2));
        }
      }
    } else { // type 0 0'1', pos=0, 1, 2
      const I ve[2] = {v0[0], v0[2]};
      bool found = m.find_simplex(1, ve, oeid);
      assert(found);

      for (const auto otriid : m.side_of(1, oeid)) {
        I tri[3];
        m.get_simplex(2, otriid, tri);

        const int pos = position12(ve, tri);
        if (pos == 0) { // 1 1'2'
          // 0 1 1'2'
          results.insert(otriid + t*n(3) + 4*m.n(3) + m.n(2));
        } else if (pos == 1) { // 0 0'2'
          // 0 0'1'2'
          results.insert(otriid + t*n(3) + 4*m.n(3) + 2*m.n(2));
        } else if (pos == 2) { // 0 0'1'
          // 0 0'1'2'
          results.insert(otriid + t*n(3) + 4*m.n(3) + 2*m.n(2));
        }
      }
    }
  }

  return results;
}

template <typename I, typename F>
int simplicial_unstructured_extruded_3d_mesh<I, F>::simplex_type(int d, I k) const
{
  const I i = mod(k, n(d)), t = std::floor((double)k / n(d));

  switch (d) {
    case 0:
      return 0;

    case 1: 
      if (i < m.n(1)) return 0;
      else if (i < 2*m.n(1)) return 1;
      else return 2;

    case 2: 
      if (i < m.n(2)) return 0;
      else if (i < 2*m.n(2)) return 1;
      else if (i < 3*m.n(2)) return 2;
      else if (i < 3*m.n(2) + m.n(1)) return 3;
      else return 4;

    case 3: 
      if (i < m.n(3)) return 0;
      else if (i < 2*m.n(3)) return 1;
      else if (i < 3*m.n(3)) return 2;
      else if (i < 4*m.n(3)) return 3;
      else if (i < 4*m.n(3) + m.n(2)) return 4;
      else if (i < 4*m.n(3) + 2*m.n(2)) return 5;
      else return 6;

    case 4:
      return i / m.n(3);

    default: 
      assert(false);
      return -1;
  }
}

}

#endif
