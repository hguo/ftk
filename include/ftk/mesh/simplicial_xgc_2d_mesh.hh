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
  void read_bfield_h5(const std::string& filename);

  I nextnode(I i) const { return nextnodes[i]; }

  const ndarray<F>& get_bfield() const { return bfield; }
  const ndarray<F>& get_psifield() const { return psifield; }

protected:
  ndarray<F> psifield, bfield;
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
  psifield(psi_),
  nextnodes(nextnodes_)
{
}

template <typename I, typename F>
void simplicial_xgc_2d_mesh<I, F>::read_bfield_h5(const std::string& filename)
{
  bfield.from_h5(filename, "/node_data[0]/values");
  bfield.set_multicomponents();
  // std::cerr << bfield.shape() << std::endl;
}

template <typename I, typename F>
std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> simplicial_xgc_2d_mesh<I, F>::from_xgc_mesh_h5(const std::string& filename)
{
  ndarray<I> triangles;
  ndarray<F> coords;
  ndarray<I> nextnodes;
  ndarray<F> psifield;

  triangles.from_h5(filename, "/cell_set[0]/node_connect_list");
  coords.from_h5(filename, "/coordinates/values");
  psifield.from_h5(filename, "/psi");
  nextnodes.from_h5(filename, "/nextnode");

  return std::shared_ptr<simplicial_xgc_2d_mesh<I, F>>(
      new simplicial_xgc_2d_mesh<I, F>(coords, triangles, psifield, nextnodes));

  // return std::shared_ptr<simplicial_unstructured_2d_mesh<I, F>>(
  //     new simplicial_unstructured_2d_mesh<I, F>(coords, triangles));
}

} // namespace ftk

#endif
