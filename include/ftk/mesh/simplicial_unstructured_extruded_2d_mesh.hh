#ifndef _FTK_simplicial_unstructured_extruded_2d_mesh_HH
#define _FTK_simplicial_unstructured_extruded_2d_mesh_HH

#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_unstructured_extruded_2d_mesh : public object {
  simplicial_unstructured_extruded_2d_mesh(std::shared_ptr<simplicial_unstructured_2d_mesh<I, F>> m_) 
    : m2(m_) {}

  virtual size_t n(int d, bool part=false) const = 0;
  size_t n_ordinal(int d, bool part=false) const { return m2->n(d, part); }
  virtual size_t n_interval(int d, bool part=false) const = 0;

  bool is_ordinal(int d, I i, bool part=false) const {
    I k = mod(i, n(d, part));
    return k < n_ordinal(d, part);
  }

  I flat_vertex_id(I i, bool part=false) const { return mod(i, m2->n(0, part)); }
  I flat_vertex_time(I i, bool part=false) const { return i / m2->n(0, part); }
  I extruded_vertex_id(I i, bool t=true, bool part=false) { return t ? i + m2->n(0, part) : i; }

  virtual int get_triangle_chi(I i) const;
  bool is_partial() const { return m2->is_partial(); }

public: // element iteration
  void element_for(int d, std::function<void(I)> f, 
      // bool part = false,
      int xl = FTK_XL_NONE, int nthreads=std::thread::hardware_concurrency(), bool affinity = false) const;

  void element_for_ordinal(int d, int t, std::function<void(I)> f,
      bool part = false,
      int xl = FTK_XL_NONE, int nthreads=std::thread::hardware_concurrency(), bool affinity = false) const;
  
  void element_for_interval(int d, int t, std::function<void(I)> f,
      bool part = false,
      int xl = FTK_XL_NONE, int nthreads=std::thread::hardware_concurrency(), bool affinity = false) const;

public: // partial base mesh
  // I lid2gid(int d, I i) const;
  // I gid2lid(int d, I i) const;

public: // mesh access
  virtual std::set<I> sides(int d, I i) const = 0;
  virtual std::set<I> side_of(int d, I i) const = 0;

  virtual I get_simplex(int d, I i, I verts[], bool part = false) const = 0;
  void get_coords(I i, F coords[]) const;

  virtual bool find_simplex(int d, I v[], I &i) const { return false; }
  virtual std::set<I> get_vertex_edge_vertex(I i) const { return {}; }

  static I mod(I val, I m);

protected:
  std::shared_ptr<simplicial_unstructured_2d_mesh<I, F>> m2;
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
void simplicial_unstructured_extruded_2d_mesh<I, F>::element_for(int d, std::function<void(I)> f, 
    int xl, int nthreads, bool affinity) const
{
  parallel_for(n(d), [&](int i) {f(i);}, xl, nthreads, affinity);
  // for (auto i = 0; i < n(d); i ++)
  //   f(i);
}

template <typename I, typename F>
void simplicial_unstructured_extruded_2d_mesh<I, F>::element_for_ordinal(
    int d, int t, std::function<void(I)> f,
    bool part,
    int xl, int nthreads, bool affinity) const
{
  parallel_for(n_ordinal(d, part), [&](int i) {
        const auto id = i + t * n(d, part);
        f(id);
      }, xl, nthreads, affinity);

  // for (auto i = 0; i < n_ordinal(d); i ++)
  //   f(i + t * n(d));
}

template <typename I, typename F>
void simplicial_unstructured_extruded_2d_mesh<I, F>::element_for_interval(int d, int t, std::function<void(I)> f,
    bool part,
    int xl, int nthreads, bool affinity) const
{
  parallel_for(n_interval(d, part), [&](int i) {
    const auto id = i + n_ordinal(d, part) + t * n(d, part);
    f(id);
  }, xl, nthreads, affinity);

  // for (auto i = 0; i < n_interval(d); i ++) {
  //   const auto id = i + n_ordinal(d) + t * n(d);
  //   f(id);
    // if (d == 2) side_of(2, id);
  // }
}

template <typename I, typename F>
void simplicial_unstructured_extruded_2d_mesh<I, F>::get_coords(I i, F coords[]) const
{
  const I k = flat_vertex_id(i),
          t = flat_vertex_time(i);
  m2->get_coords(k, coords);

  coords[3] = t;
  // coords[ m2->ncoords() ] = t; // last coordinate
}

template <typename I, typename F>
int simplicial_unstructured_extruded_2d_mesh<I, F>::get_triangle_chi(I i) const
{
  const int d = 2;
  const I k = mod(i, n(d)); // t = std::floor(double(k) / n(d));
  if (k < n_ordinal(d)) {
    return m2->get_triangle_chi(k);
  } else 
    return 0; // no chirality
}

}

#endif
