#ifndef _HYPERMESH_simplicial_unstructured_extruded_2d_mesh_implicit_HH
#define _HYPERMESH_simplicial_unstructured_extruded_2d_mesh_implicit_HH

#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_unstructured_extruded_2d_mesh_implicit : public simplicial_unstructured_extruded_2d_mesh<I, F> {
  simplicial_unstructured_extruded_2d_mesh_implicit(std::shared_ptr<simplicial_unstructured_2d_mesh<I, F>> m) 
    : simplicial_unstructured_extruded_2d_mesh<I, F>(m) {}

  size_t n(int d, bool part=false) const;
  size_t n_interval(int d, bool part=false) const;

  int face_type(I i) const; // legacy
  int simplex_type(int d, I i, bool part=false) const;

public: // mesh access
  std::set<I> sides(int d, I i) const;
  std::set<I> side_of(int d, I i) const;

  I get_simplex(int d, I i, I verts[], bool part = false) const;

  bool find_simplex(int d, I v[], I &i) const;
 
  std::set<I> get_vertex_edge_vertex(I i) const;

protected: // legacy code to produce an explicit extrusion
  void extrude();

  // ftk::ndarray<I> tets, // {4, n(3)}
  //                 tris, // {3, n(2)}
  //                 edges;// {2, n(1)}
};

/////

template <typename I, typename F>
void simplicial_unstructured_extruded_2d_mesh_implicit<I, F>::extrude() 
{
  fprintf(stderr, "3d mesh initialized, #tet=%zu, #tri=%zu, #edge=%zu, #vert=%zu\n",
      n(3), n(2), n(1), n(0));

#if 0
  tets.reshape({4, n(3)});
  
  // tets
  for (auto i = 0; i < this->m2->n(2); i ++) {
    I tri[3];
    this->m2->get_triangle(i, tri);
    std::sort(tri, tri+3);

    // tets
    tets(0, i) = tri[0]; // type I: bottom tet
    tets(1, i) = tri[1];
    tets(2, i) = tri[2];
    tets(3, i) = tri[2] + this->m2->n(0);

    tets(0, i+this->m2->n(2)) = tri[0]; // type II: center tet
    tets(1, i+this->m2->n(2)) = tri[1];
    tets(2, i+this->m2->n(2)) = tri[1] + this->m2->n(0);
    tets(3, i+this->m2->n(2)) = tri[2] + this->m2->n(0);

    tets(0, i+2*this->m2->n(2)) = tri[0]; // type III: top tet
    tets(1, i+2*this->m2->n(2)) = tri[0] + this->m2->n(0);
    tets(2, i+2*this->m2->n(2)) = tri[1] + this->m2->n(0);
    tets(3, i+2*this->m2->n(2)) = tri[2] + this->m2->n(0);
  }

  // unique triangles
  tris.reshape({3, n(2)}); 
  for (auto i = 0; i < this->m2->n(2); i ++) { // per prism
    I tri[3];
    this->m2->get_triangle(i, tri);
    
    tris(0, i) = tri[0]; // type I: base triangle
    tris(1, i) = tri[1];
    tris(2, i) = tri[2];

    std::sort(tri, tri+3); // the following extrusions need sorted vertices
    tris(0, i + this->m2->n(2)) = tri[0]; // type II: prism lower triangle
    tris(1, i + this->m2->n(2)) = tri[1];
    tris(2, i + this->m2->n(2)) = tri[2] + this->m2->n(0);
    
    tris(0, i + this->m2->n(2)*2) = tri[0]; // type III: prism upper triangle
    tris(1, i + this->m2->n(2)*2) = tri[1] + this->m2->n(0);
    tris(2, i + this->m2->n(2)*2) = tri[2] + this->m2->n(0);
  }
  for (auto i = 0; i < this->m2->n(1); i ++) { // per edge 
    I v[2];
    this->m2->get_edge(i, v);
    if (v[0] > v[1]) std::swap(v[0], v[1]); // sorted

    tris(0, i + this->m2->n(2)*3) = v[0]; // type IV: edge lower triangle
    tris(1, i + this->m2->n(2)*3) = v[1];
    tris(2, i + this->m2->n(2)*3) = v[1] + this->m2->n(0);
    
    tris(0, i + this->m2->n(2)*3 + this->m2->n(1)) = v[0]; // type V: edge upper triangle
    tris(1, i + this->m2->n(2)*3 + this->m2->n(1)) = v[0] + this->m2->n(0);
    tris(2, i + this->m2->n(2)*3 + this->m2->n(1)) = v[1] + this->m2->n(0);
  }

  // unique edges
  edges.reshape({2, n(1)});
  for (auto i = 0; i < this->m2->n(1); i ++) {
    I v[2];
    this->m2->get_edge(i, v);
    if (v[0] > v[1]) std::swap(v[0], v[1]);

    edges(0, i) = v[0];
    edges(1, i) = v[1];
    
    edges(0, i+this->m2->n(1)) = v[0];
    edges(1, i+this->m2->n(1)) = v[1] + this->m2->n(0);
  }
  for (auto i = 0; i < this->m2->n(0); i ++) {
    edges(0, i+this->m2->n(1)*2) = i; 
    edges(1, i+this->m2->n(1)*2) = i + this->m2->n(0);
  }
#endif
}

template <typename I, typename F>
size_t simplicial_unstructured_extruded_2d_mesh_implicit<I, F>::n(int d, bool part) const
{
  if (d == 0) {
    return this->m2->n(0, part); // unique vertices in each interval
  } else if (d == 1) {
    return this->m2->n(1, part) * 2 + this->m2->n(0, part); // unique edges in each interval
  } else if (d == 2) {
    return this->m2->n(2, part) // the bottom triangle of the prism
      + this->m2->n(2, part) * 2 // two triangles inside prisms
      + this->m2->n(1, part) * 2;// unique triangles extruded from each edge
  } else if (d == 3) { 
    return 3 * this->m2->n(2, part); // 3 tets per prism
  } else return 0;
}

template <typename I, typename F>
size_t simplicial_unstructured_extruded_2d_mesh_implicit<I, F>::n_interval(int d, bool part) const
{
  if (d == 1) return this->m2->n(1, part) + this->m2->n(0, part); // return this->m2->n(1)*2;
  else if (d == 2) return this->m2->n(2, part)*2 + this->m2->n(1, part)*2;
  else if (d == 3) return this->m2->n(2, part)*3;
  else return 0;
}

#if 0
template <typename I, typename F>
I simplicial_unstructured_extruded_2d_mesh_implicit<I, F>::lid2gid(int d, I lid) const
{
  const bool part = true;
  const I i = mod(k, n(d, part)), t = std::floor(double(k) / n(d, part));
  // const I offset = t * n(0, part);
  const int type = simplex_type(d, i, part);

  if (d == 0) {
    return this->m2->lid2gid(d, i) + t * n(d);
  } else if (d == 1) { // derive flat vert/edge lid and translate to gid
    if (i < this->m2->n(1, part)) { // transloate local vertex id

  } else if (d == 2) { // derive flat edge/tri lid and translate to gid

  } else if (d == 3) { // derive flat tet lid and translate to gid

  } else 
    return -1; // assert false
}
#endif

template <typename I, typename F>
I simplicial_unstructured_extruded_2d_mesh_implicit<I, F>::get_simplex(int d, I k, I verts[], bool part) const
{
  const I i = this->mod(k, n(d, part)), t = std::floor(double(k) / n(d, part));
  const I offset = t * n(0); // 0d gid offset for verts
  const int type = simplex_type(d, i, part);
  I gid;

  if (d == 0) {
    if (part) verts[0] = this->m2->lid2gid(0, i);
    else verts[0] = i; // + offset;
    gid = verts[0];
  } else if (d == 1) {
    if (type < 2) {
      I edge[2];
      I eid = this->m2->get_simplex(1, i % this->m2->n(1, part), edge, part);

      switch (type) {
        case 0: // 0 1
          verts[0] = edge[0];
          verts[1] = edge[1];
          gid = eid;
          break;

        case 1: // 0 1'
          verts[0] = edge[0];
          verts[1] = edge[1] + this->m2->n(0);
          gid = eid + this->m2->n(1);
      }
    } else {
      const I vid = part ? this->m2->lid2gid(0, i) : i;
      verts[0] = vid - 2*this->m2->n(1);
      verts[1] = vid - 2*this->m2->n(1) + this->m2->n(0);
      gid = vid; // + 2 * this->m2->n(1);
    }
  } else if (d == 2) {
    if (type < 3) {
      I tri[3];
      I tid = this->m2->get_simplex(2, i % this->m2->n(2, part), tri, part);
      switch (type) {
        case 0: // 0 1 2
          verts[0] = tri[0];
          verts[1] = tri[1];
          verts[2] = tri[2];
          gid = tid;
          break;

        case 1: // 0 1 2'
          verts[0] = tri[0];
          verts[1] = tri[1];
          verts[2] = tri[2] + this->m2->n(0);
          gid = tid + this->m2->n(2);
          break;

        case 2: // 0 1'2'
          verts[0] = tri[0];
          verts[1] = tri[1] + this->m2->n(0);
          verts[2] = tri[2] + this->m2->n(0);
          gid = tid + 2*this->m2->n(2);
          break;

        default: assert(false);
      }
    } else {
      I edge[2];
      I eid = this->m2->get_simplex(1, (i - 3*this->m2->n(2, part)) % this->m2->n(1, part), edge, part);
      switch (type) {
        case 3: 
          verts[0] = edge[0];
          verts[1] = edge[1];
          verts[2] = edge[1] + this->m2->n(0);
          gid = eid + 3*this->m2->n(2);
          break;

        case 4:
          verts[0] = edge[0];
          verts[1] = edge[0] + this->m2->n(0);
          verts[2] = edge[1] + this->m2->n(0);
          gid = eid + 3*this->m2->n(2) + this->m2->n(1);
          break;

        default: assert(false);
      }
    }
  } else if (d == 3) {
    I tri[3];
    I tid = this->m2->get_simplex(2, i % this->m2->n(2, part), tri, part);

    switch (type) {
      case 0: // type I: 0 1 2 2'
        verts[0] = tri[0];
        verts[1] = tri[1];
        verts[2] = tri[2];
        verts[3] = tri[2] + this->m2->n(0);
        gid = tid;
        break;

      case 1: // type II: 0 1 1'2'
        verts[0] = tri[0];
        verts[1] = tri[1];
        verts[2] = tri[1] + this->m2->n(0);
        verts[3] = tri[2] + this->m2->n(0);
        gid = tid + this->m2->n(2);
        break;
      
      case 2: // type III: 0 0'1'2'
        verts[0] = tri[0];
        verts[1] = tri[0] + this->m2->n(0);
        verts[2] = tri[1] + this->m2->n(0);
        verts[3] = tri[2] + this->m2->n(0);
        gid = tid + 2*this->m2->n(2);
        break;

      default:
        assert(false);
    }
  }
  
  for (int i = 0; i < d+1; i ++)
    verts[i] += offset;

  gid = gid + t * n(d);
  return gid;

#if 0
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
#endif
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_extruded_2d_mesh_implicit<I, F>::sides(int d, I k) const 
{
  const I i = this->mod(k, n(d)), t = std::floor((double)k / n(d));
  std::set<I> results;

  I v[4];
  get_simplex(3, i, v);

  if (d == 3) { // currently only support d==3
    const int type = i / this->m2->n(2);
    // fprintf(stderr, "tet_type=%d\n", type);
    if (type == 0) { // bottom 0-1-2-2'
      // 0-1-2, type I
      // 0-1-2', type II
      // 1-2-2', type IV
      // 0-2-2', type IV
      {
        const I ot[3] = {v[0], v[1], v[2]}; // triangle in the orignal mesh
        I otid; // triangle id in the original mesh
        bool found = this->m2->find_triangle(ot, otid);
        assert(found);

        results.insert( otid + t*n(2) );
        results.insert( otid + t*n(2) + this->m2->n(2) );
      }
      {
        const I oe[2] = {v[1], v[2]};
        I oeid; // edge id in the orignal mesh
        bool found = this->m2->find_edge(oe, oeid);
        assert(found);
        results.insert( oeid + this->m2->n(2)*3 + t*n(2));
      }
      {
        const I oe[2] = {v[0], v[2]};
        I oeid; // edge id in the orignal mesh
        bool found = this->m2->find_edge(oe, oeid);
        assert(found);
        results.insert( oeid + this->m2->n(2)*3 + t*n(2));
      }
    } else if (type == 1) { // center 0-1-1'-2'
      // 0-1-2', type II
      // 0-1'-2', type III
      // 0-1-1', type IV
      // 1-1'-2', type V
      {
        const I ot[3] = {v[0], v[1], this->mod(v[3], this->m2->n(0))}; // triangle in the orignal mesh
        I otid; // triangle id in the original mesh
        bool found = this->m2->find_triangle(ot, otid);
        assert(found);

        results.insert( otid + this->m2->n(2) + t*n(2) );
        results.insert( otid + this->m2->n(2)*2 + t*n(2) );
      }
      {
        const I oe[2] = {v[0], v[1]};
        I oeid; // edge id in the orignal mesh
        bool found = this->m2->find_edge(oe, oeid);
        assert(found);
        results.insert( oeid + this->m2->n(2)*3 + t*n(2) );
      }
      {
        const I oe[2] = {v[1], this->mod(v[3], n(0))}; 
        I oeid; // edge id in the orignal mesh
        bool found = this->m2->find_edge(oe, oeid);
        assert(found);
        results.insert( oeid + this->m2->n(2)*3 + this->m2->n(1) + t*n(2) );
      }
    } else if (type == 2) { // top 0-0'-1'-2'
      // 0-1-2, t+1, type I
      // 0-1'-2', type III
      // 0-0'-1', type V
      // 0-0'-2', type V
      {
        const I ot[3] = {v[0], this->mod(v[2], this->m2->n(0)), this->mod(v[3], this->m2->n(0))}; // triangle in the orignal mesh
        I otid; // triangle id in the original mesh
        bool found = this->m2->find_triangle(ot, otid);
        assert(found);

        results.insert( otid + (t+1)*n(2) );
        results.insert( otid + this->m2->n(2)*2 + t*n(2) );
      }
      {
        const I oe[2] = {v[0], this->mod(v[2], n(0))}; 
        I oeid; // edge id in the orignal mesh
        bool found = this->m2->find_edge(oe, oeid);
        assert(found);
        results.insert( oeid + this->m2->n(2)*3 + this->m2->n(1) + t*n(2) );
      }
      {
        const I oe[2] = {v[0], this->mod(v[3], this->m2->n(0))}; 
        I oeid; // edge id in the orignal mesh
        bool found = this->m2->find_edge(oe, oeid);
        assert(found);
        results.insert( oeid + this->m2->n(2)*3 + this->m2->n(1) + t*n(2) );
      }
    }
  }

  return results;
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_extruded_2d_mesh_implicit<I, F>::side_of(int d, I k) const 
{
  const I i = this->mod(k, n(d)), t = std::floor((double)k / n(d));
  // fprintf(stderr, "k=%d, i=%d, t=%d, n(d)=%d\n", k, i, t, n(d));
  std::set<I> results;

  if (d == 2) { // currently only support d==2
    int v[3];
    get_simplex(2, k, v);

    // fprintf(stderr, "k=%d, i=%d, t=%d, v=%d, %d, %d\n", k, i, t, v[0], v[1], v[2]);
    // fprintf(stderr, "face=%d,%d type=%d\n", k, i, face_type(k));

    // determine triangle type by id
    const I v0[3] = {this->mod(v[0], this->m2->n(0)), this->mod(v[1], this->m2->n(0)), this->mod(v[2], this->m2->n(0))};
    if (i < 3*this->m2->n(2)) {
      I otid; // triangle id in the original mesh
      bool found = this->m2->find_triangle(v0, otid); 
      assert(found);
      
      if (i < this->m2->n(2)) { // type I: bottom tet
        // t-1: 0-0'-1'-2' (type III tet)
        // t:   0-1-2-2' (type I tet)
        results.insert( otid + 2*this->m2->n(2) + (t-1)*n(3) );
        results.insert( otid            + t    *n(3) );
      } else if (i < 2*this->m2->n(2)) { // type II: prism lower triangle
        // t:   0-1-2-2' (type I tet)
        // t:   0-1-1'-2' (type II tet)
        results.insert( otid          + t*n(3) );
        results.insert( otid + this->m2->n(2) + t*n(3) );
      } else { // if (i < 3*this->m2->n(2)) { // type III: prism upper triangle
        // t:   0-1-1'-2'  (type II tet)
        // t:   0-0'-1'-2' (type III tet)
        results.insert( otid +   this->m2->n(2) + t*n(3) );
        results.insert( otid + 2*this->m2->n(2) + t*n(3) );
      }
    } else {
      I oeid; // edge id in the original mesh
      I ve[2] = {v0[0], v0[2]};
      bool found = this->m2->find_edge(ve, oeid);
      assert(found);

      const bool lower = i < 3*this->m2->n(2) + this->m2->n(1);

      const auto triangles = this->m2->side_of(1, oeid);
      for (const auto triangle : triangles) {
        I vc[3];
        this->m2->get_triangle(triangle, vc);

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
            results.insert(triangle + this->m2->n(2) + t*n(3)); // center tet
        } else { // type V: edge upper triangle
          if (pos == 0) 
            results.insert(triangle + this->m2->n(2) + t*n(3)); // center tet
          else if (pos == 1)
            results.insert(triangle + 2*this->m2->n(2) + t*n(3)); // top tet
          else if (pos == 2)
            results.insert(triangle + 2*this->m2->n(2) + t*n(3)); // top tet
        }
      }

#if 0
      int tids[2] = {-1}; // triangle ids in array
      int j = 0;
      for (auto tri : triangles) {
        tids[j++] = tri;

        I ve1[2];
        this->m2->get_edge(oeid, ve1);

        int vv[3];
        this->m2->get_triangle(tri, vv);
        fprintf(stderr, "ve=%d, %d, oeid=%d, ve1=%d, %d, vv=%d, %d, %d, tri=%d\n", 
            ve[0], ve[1], oeid, 
            ve1[0], ve1[1], 
            vv[0], vv[1], vv[2], tri);
      }
#endif

#if 0
      // find triangle v0-v1-v2, denote the other triangle as vn-v0-v1
      if (i < 3*this->m2->n(2) + this->m2->n(1)) { // type IV: edge lower triangle
        // type II tet for the large triangle
        // type I tet for the small triangle
        if (large_triangle >= 0) 
          results.insert( large_triangle + this->m2->n(2) + t*n(3) );
        if (small_triangle >= 0) 
          results.insert( small_triangle          + t*n(3) );
      } else { // type V: edge upper triangle
        // type III tet for the large triangle
        // type II tet for the small triangle
        if (large_triangle >= 0)
          results.insert( large_triangle + 2*this->m2->n(2) + t*n(3) );
        if (small_triangle >= 0) 
          results.insert( small_triangle +   this->m2->n(2) + t*n(3) );
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
int simplicial_unstructured_extruded_2d_mesh_implicit<I, F>::face_type(I k) const
{
  return simplex_type(2, k);
#if 0
  const I i = mod(k, n(2));

  if (i < this->m2->n(2)) return 0;  // I
  else if (i < 2*this->m2->n(2)) return 1; // II
  else if (i < 3*this->m2->n(2)) return 2; // III 
  else if (i < 3*this->m2->n(2) + this->m2->n(1)) return 3; // IV
  else if (i < 3*this->m2->n(2) + 2*this->m2->n(1)) return 4; // V
  else return -1;
#endif
}

template <typename I, typename F>
int simplicial_unstructured_extruded_2d_mesh_implicit<I, F>::simplex_type(int d, I k, bool part) const
{
  const I i = this->mod(k, n(d, part)), t = std::floor((double)k / n(d, part));

  switch (d) {
    case 0:
      return 0;

    case 1: 
      if (i < this->m2->n(1, part)) return 0;
      else if (i < 2*this->m2->n(1, part)) return 1;
      else return 2;

    case 2: 
      if (i < this->m2->n(2, part)) return 0;
      else if (i < 2*this->m2->n(2, part)) return 1;
      else if (i < 3*this->m2->n(2, part)) return 2;
      else if (i < 3*this->m2->n(2, part) + this->m2->n(1, part)) return 3;
      else return 4;

    case 3:
      return i / this->m2->n(2, part);

    default: 
      assert(false);
      return -1;
  }
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_extruded_2d_mesh_implicit<I, F>::get_vertex_edge_vertex(I i) const
{
  const I k = this->flat_vertex_id(i),
          t = this->flat_vertex_time(i);

  std::set<I> verts0 = this->m2->get_vertex_edge_vertex(k);
  std::set<I> verts;
  for (const auto v0 : verts0)
    verts.insert( t * this->m2->n(0) + v0 ); // neighbors in the same timestep
  verts.insert( (t+1) * this->m2->n(0) + k );
  verts.insert( (t-1) * this->m2->n(0) + k );

  return verts;
}

template <typename I, typename F>
bool simplicial_unstructured_extruded_2d_mesh_implicit<I, F>::find_simplex(int d, I v[], I &id) const
{
  I v0[d+1];
  int t[d+1], t0[d+1], mint = std::numeric_limits<int>::max();
  std::set<I> vset;

  for (int i = 0; i < d+1; i ++) {
    v0[i] = this->flat_vertex_id(v[i]);
    t[i] = this->flat_vertex_time(v[i]);
    mint = std::min(mint, t[i]);
    vset.insert(v0[i]);
  }
  // fprintf(stderr, "***** mint=%d, t=%d, %d, %d, %d, v0=%d, %d, %d, %d\n", 
  //     mint, t[0], 
  //     t[1], t[2], t[3], 
  //     v0[0], v0[1], v0[2], v0[3]);

  for (int i = 0; i < d+1; i ++)
    t0[i] = t[i] - mint;
    
  int w[4];
  int k = 0;
  for (auto v : vset) 
    w[k++] = v;

  if (d == 3) {
    if (!this->m2->find_triangle(w, id)) return false;
    if (t0[1] == 0 && t0[2] == 0) // type I
    {} // nothing to do
    else if (t0[1] == 0 && t0[2] == 1) // type II
      id += this->m2->n(2);
    else if (t0[1] == 1 && t0[2] == 1) // type III
      id += 2*this->m2->n(2);
    else {
      fprintf(stderr, "FAAAATAL: v=%d, %d, %d, %d, v0=%d, %d, %d, %d, t0=%d, %d, %d, %d\n", 
          v[0], v[1], v[2], v[3], 
          v0[0], v0[1], v0[2], v0[3],
          t0[0], t0[1], t0[2], t0[3]);
      assert(false);
    }
    id += mint * n(3);
    return true;
  } else if (d == 2) {
    if (vset.size() == 3) {
      if (!this->m2->find_triangle(w, id)) return false;
      if (t0[1] == 0 && t0[2] == 0) // type I
      {} // nothing to do
      else if (t0[1] == 0 && t0[2] == 1) // type II
        id += this->m2->n(2);
      else if (t0[1] == 1 && t0[2] == 1) // type III
        id += 2*this->m2->n(2);
      else assert(false);
    } else if (vset.size() == 2) {
      if (!this->m2->find_edge(w, id)) return false;
      if (t0[1] == 0) // type IV
        id += 3*this->m2->n(2);
      else if (t0[1] == 1) // type V
        id += 3*this->m2->n(2) + this->m2->n(1);
      else assert(false);
    } else assert(false);
    id += mint * n(2);
    return true;
  } else if (d == 1) {
    if (vset.size() == 2) {
      if (!this->m2->find_edge(w, id)) return false;
      // fprintf(stderr, "***** t0=%d, %d, w=%d, %d, id=%d\n", t0[0], t0[1], w[0], w[1], id);

      if (t0[1] == 0) // type I
      {} // nothing to do
      else if (t0[1] == 1) // type II
        id += this->m2->n(1);
      else assert(false);
    } else if (vset.size() == 1) 
      id = *vset.begin() + this->m2->n(1)*2;
    else assert(false);
    id += mint * n(1);
    return true;
  }
  return false;
}

}

#endif
