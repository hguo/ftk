#ifndef _FTK_XGC_2D_MESH_HH
#define _FTK_XGC_2D_MESH_HH

#include <ftk/ftk_config.hh>
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>
#include <ftk/mesh/point_locator_2d_quad.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/print.hh>
#include <ftk/io/xgc_units.hh>

namespace ftk {

template <typename I=int, typename F=double>
struct simplicial_xgc_2d_mesh : public simplicial_unstructured_2d_mesh<I, F> {
  simplicial_xgc_2d_mesh(
      const ndarray<F>& coords, 
      const ndarray<I>& triangles,
      const ndarray<F>& psi,
      const ndarray<I>& nextnodes);

  void initialize_point_locator();

  static std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> from_xgc_mesh_h5(const std::string& filename);
  void read_bfield_h5(const std::string& filename);
  bool read_units_m(const std::string& filename) { return units.read(filename); }

  I nextnode(I i) const { return nextnodes[i]; }

  const xgc_units_t& get_units() const { return units; }
  const ndarray<I>& get_nextnodes() const { return nextnodes; }
  const ndarray<F>& get_bfield() const { return bfield; }
  const ndarray<F>& get_psifield() const { return psifield; }

  bool eval_b(const F x[], F b[]) const; // get magnetic field value at x
  bool eval_f(const F rzp[3], F f[2]) const;
  F eval_psi(const F x[]) const; // get psi value at x

  void magnetic_map(F rzp[3], F phi_end, int nsteps=100) const;

protected:
  xgc_units_t units;
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
void simplicial_xgc_2d_mesh<I, F>::magnetic_map(F rzp[3], F phi_end, int nsteps) const
{
  const F delta = (phi_end - rzp[2]) / nsteps;

  // rk1
  for (int k = 0; k < nsteps; k ++) {
    F B[3];
    bool succ = eval_b(rzp, B);

    if (succ) {
      const F r = rzp[0], z = rzp[1], phi = rzp[2];
      rzp[0] += delta * r * B[0] / B[2];
      rzp[1] += delta * r * B[1] / B[2];
    }
    rzp[2] += delta;
    // fprintf(stderr, "integrating %f, %f, %f\n", rzp[0], rzp[1], rzp[2]);
  }
}

template <typename I, typename F>
void simplicial_xgc_2d_mesh<I, F>::initialize_point_locator()
{
  this->locator.reset( new point_locator_2d_quad<I, F>(*this) );
}

template <typename I, typename F>
bool simplicial_xgc_2d_mesh<I, F>::eval_b(const F x[], F b[]) const
{
  F mu[3];
  I tid = this->locator->locate(x, mu);
  if (tid < 0) {
    b[0] = b[1] = b[2] = F(0); 
    return false;
  } else {
    F B[3][3];
    I tri[3];
    this->get_simplex(2, tid, tri);
    for (int i = 0; i < 3; i ++) 
      for (int j = 0; j < 3; j ++) 
        B[i][j] = bfield(j, tri[i]);
    lerp_s2v3(B, mu, b);
    
    // print3x3("B", B);
    // fprintf(stderr, "mu=%f, %f, %f, b=%f, %f, %f\n", 
    //     mu[0], mu[1], mu[2], b[0], b[1], b[2]);
    return true;
  }
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
