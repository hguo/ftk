#ifndef _FTK_XGC_2D_MESH_HH
#define _FTK_XGC_2D_MESH_HH

#include <ftk/config.hh>
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
  void initialize_roi(double psin_min = 0.9, double psin_max = 1.05, 
      double theta_min_deg = -60.0, double theta_max_deg = 60.0);

  static std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> from_xgc_mesh_h5(const std::string& filename);
  static std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> from_xgc_mesh_adios2(diy::mpi::communicator comm, const std::string& filename);
  void read_bfield_h5(const std::string& filename);
  bool read_units_m(const std::string& filename) { return units.read(filename); }

  I nextnode(I i) const { return nextnodes[i]; }
  F psin(I i) const { return psifield[i] / units.psi_x; } // get normalized psi

  const xgc_units_t& get_units() const { return units; }
  const ndarray<I>& get_nextnodes() const { return nextnodes; }
  const ndarray<F>& get_bfield() const { return bfield; }
  const ndarray<F>& get_psifield() const { return psifield; }
  const std::vector<I>& get_roi_nodes() const { return roi_nodes; }

  bool eval_b(const F x[], F b[]) const; // get magnetic field value at x
  // bool eval_f(const F rzp[3], F f[2]) const;
  F eval_psi(const F x[]) const; // get psi value at x

  void magnetic_map(F rzp[3], F phi_end, int nsteps=100) const;

protected:
  xgc_units_t units;
  ndarray<F> psifield, bfield;
  ndarray<I> nextnodes;
  ndarray<I> roi;

  std::vector<I> roi_nodes;
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
  bfield.read_h5(filename, "/node_data[0]/values");
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
std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> simplicial_xgc_2d_mesh<I, F>::from_xgc_mesh_adios2(
    diy::mpi::communicator comm, 
    const std::string& filename)
{
  ndarray<I> triangles;
  ndarray<F> coords;
  ndarray<I> nextnodes;
  ndarray<F> psifield;

  triangles.from_bp(comm, filename, "/cell_set[0]/node_connect_list");
  coords.from_bp(comm, filename, "/coordinates/values");
  psifield.from_bp(comm, filename, "psi");
  nextnodes.from_bp(comm, filename, "nextnode");

  triangles.reshape(triangles.dim(1), triangles.dim(0));
  coords.reshape(coords.dim(1), coords.dim(0));

  std::cerr << triangles.shape() << std::endl;
  std::cerr << coords.shape() << std::endl;
  std::cerr << nextnodes.shape() << std::endl;
  std::cerr << psifield.shape() << std::endl;

  return std::shared_ptr<simplicial_xgc_2d_mesh<I, F>>(
      new simplicial_xgc_2d_mesh<I, F>(coords, triangles, psifield, nextnodes));
}

template <typename I, typename F>
std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> simplicial_xgc_2d_mesh<I, F>::from_xgc_mesh_h5(const std::string& filename)
{
  ndarray<I> triangles;
  ndarray<F> coords;
  ndarray<I> nextnodes;
  ndarray<F> psifield;

  triangles.read_h5(filename, "/cell_set[0]/node_connect_list");
  coords.read_h5(filename, "/coordinates/values");
  psifield.read_h5(filename, "/psi");
  nextnodes.read_h5(filename, "/nextnode");

  return std::shared_ptr<simplicial_xgc_2d_mesh<I, F>>(
      new simplicial_xgc_2d_mesh<I, F>(coords, triangles, psifield, nextnodes));

  // return std::shared_ptr<simplicial_unstructured_2d_mesh<I, F>>(
  //     new simplicial_unstructured_2d_mesh<I, F>(coords, triangles));
}

template <typename I, typename F>
void simplicial_xgc_2d_mesh<I, F>::initialize_roi(
    double psin_min, double psin_max, double theta_min_deg, double theta_max_deg)
{
  const F r0 = units.eq_axis_r, 
          z0 = units.eq_axis_z;
  const F psi_x = units.psi_x;
  const F theta_min_rad = theta_min_deg * M_PI / 180, 
          theta_max_rad = theta_max_deg * M_PI / 180;

  roi.reshape(nextnodes.shape());
  roi_nodes.clear();

  for (int i = 0; i < psifield.size(); i ++) {
    const F r = this->vertex_coords(0, i), 
            z = this->vertex_coords(1, i);
    
    const F theta = atan2(z - z0, r - r0);
    if (theta < theta_min_rad || theta > theta_max_rad) continue;

    const F psin = psifield[i] / psi_x;
    if (psin < psin_min || psin > psin_max) continue;

    roi_nodes.push_back(i);
    roi[i] = 1;
  }

  // fprintf(stderr, "#xgc_roi_nodes=%zu\n", roi_nodes.size());
  // this->array_to_vtu("xgc-roi.vtu", "roi", roi);
}

} // namespace ftk

#endif
