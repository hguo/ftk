#ifndef _FTK_XGC_2D_MESH_HH
#define _FTK_XGC_2D_MESH_HH

#include <ftk/ftk_config.hh>
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_xgc_2d_mesh : public simplicial_unstructured_2d_mesh<I, F> {
  simplicial_xgc_2d_mesh(
      const ndarray<F>& coords, 
      const ndarray<I>& triangles,
      const ndarray<F>& psi,
      const ndarray<I>& nextnodes);

  static std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> from_xgc_mesh_h5(const std::string& filename);

  I nextnode(I i) const { return nextnodes[i]; }

public:
#if FTK_HAVE_VTK
#endif

protected:
  ndarray<F> psi;
  ndarray<I> nextnodes;
};
/////////
  
template <typename I, typename F>
simplicial_xgc_2d_mesh<I, F>::simplicial_xgc_2d_mesh(
    const ndarray<F>& coords, 
    const ndarray<I>& triangles,
    const ndarray<F>& psi_,
    const ndarray<I>& nextnodes_) : 
  simplicial_unstructured_2d_mesh<I, F>(coords, triangles), 
  psi(psi_),
  nextnodes(nextnodes_)
{
}

template <typename I, typename F>
std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> simplicial_xgc_2d_mesh<I, F>::from_xgc_mesh_h5(const std::string& filename)
{
  ndarray<I> triangles;
  ndarray<F> coords;
  ndarray<I> nextnodes;
  ndarray<F> psi;

  triangles.from_h5(filename, "/cell_set[0]/node_connect_list");
  coords.from_h5(filename, "/coordinates/values");
  psi.from_h5(filename, "/psi");
  nextnodes.from_h5(filename, "/nextnode");

  return std::shared_ptr<simplicial_xgc_2d_mesh<I, F>>(
      new simplicial_xgc_2d_mesh<I, F>(coords, triangles, psi, nextnodes));

  // return std::shared_ptr<simplicial_unstructured_2d_mesh<I, F>>(
  //     new simplicial_unstructured_2d_mesh<I, F>(coords, triangles));
}

} // namespace ftk

#endif
