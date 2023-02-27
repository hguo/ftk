#ifndef _FTK_MESH_UNSTRUCTURED_PERIODIC_2D_HH
#define _FTK_MESH_UNSTRUCTURED_PERIODIC_2D_HH

#include <ftk/config.hh>
#include <ftk/mesh/simplicial_unstructured_mesh.hh>
#include <ftk/utils/string.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_unstructured_periodic_2d_mesh : public object {
  simplicial_unstructured_periodic_2d_mesh(std::shared_ptr<simplicial_unstructured_3d_mesh<I, F> m_) : m(m_) {}

  size_t n(int d) const;

protected:
  // the base mesh with identical lower and upper 2D meshes
  std::shared_ptr<simplicial_unstructured_3d_mesh<I, F>> m3; 
  std::shared_ptr<simplicial_unstructured_2d_mesh<I, F>> m2; 
};

/////////
template <typename I, typename F> 
size_t simplicial_unstructured_periodic_2d_mesh<I, F>::n(int d) const // number of elements per unit period
{
  return m3->n(d) - m2->n(d);
}

template <typename I, typename F>
I simplicial_unstructured_periodic_2d_mesh<I, F>::get_simplex(int d, I k, I verts[]) const
{
  const I i = mod(k, n(d, part)), t = std::floor(double(k) / n(d, part));
  const I offset = t * n(0); // 0d gid offset for verts

  m3->get_simplex(d, k, verts);
  for (int i = 0; i < d; i ++)
    verts[i] += offset;

  I gid = gid + t * n(d);
  return gid;
}

template <typename I, typename F>
void simplicial_unstructured_periodic_2d_mesh<I, F>::get_coords(I i, F coords[]) const
{
  const I k = flat_vertex_id(i),
          t = flat_vertex_time(i);
  m.get_coords(k, coords);

  coords[3] = t;
}

template <typename I, typename F>
std::set<I> simplicial_unstructured_periodic_2d_mesh<I, F>::sides(int d, I k) const 
{
  const I i = mod(k, n(d)), t = std::floor((double)k / n(d));
  std::set<I> results;

  std::set<I> m3sides = m3->sides(d, i);

  // TODO: not implemented yet

  return results;
}

#endif
