#ifndef _FTK_XGC_2D_MESH_HH
#define _FTK_XGC_2D_MESH_HH

#include <ftk/config.hh>
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>
#include <ftk/mesh/point_locator_2d_quad.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/fmod.hh>
#include <ftk/numeric/rk4.hh>
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
  
  simplicial_xgc_2d_mesh(
      const ndarray<F>& coords, 
      const ndarray<I>& triangles) : simplicial_unstructured_2d_mesh<I, F>(coords, triangles) {}

  void initialize_point_locator();
  void initialize_roi(double psin_min = 0.9, double psin_max = 1.05, 
      double theta_min_deg = -60.0, double theta_max_deg = 60.0);

  static std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> from_xgc_mesh_file(const std::string& filename, diy::mpi::communicator comm = MPI_COMM_WORLD);
  static std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> from_xgc_mesh_h5(const std::string& filename);
  static std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> from_xgc_mesh_bp(const std::string& filename, diy::mpi::communicator comm = MPI_COMM_WORLD);

  void read_oneddiag(const std::string& filename, diy::mpi::communicator comm = MPI_COMM_WORLD);
  void read_bfield(const std::string& filename, diy::mpi::communicator comm = MPI_COMM_WORLD);
  bool read_units_m(const std::string& filename) { return units.read(filename); }

  I nextnode(I i) const { return nextnodes[i]; }
  F psin(I i) const { return psifield[i] / units.psi_x; } // get normalized psi

  void set_units(const xgc_units_t& u) { units = u; }
  const xgc_units_t& get_units() const { return units; }
  const ndarray<I>& get_nextnodes() const { return nextnodes; }
  const ndarray<F>& get_bfield() const { return bfield; }
  const ndarray<F>& get_bfield0() const { return bfield0; }
  const ndarray<F>& get_curl_bfield0() const { return curl_bfield0; }
  const ndarray<F>& get_psifield() const { return psifield; }
  ndarray<F> get_psinfield() const { return psifield * (1.0 / units.psi_x); }
  const std::vector<I>& get_roi_nodes() const { return roi_nodes; }

  bool eval_b(const F x[], F b[]) const; // get magnetic field value at x
  bool eval_total_B(const ndarray<F>& totalB, const F rzp[3], F b[3]) const;
  // bool eval_f(const F rzp[3], F f[2]) const;
  F eval_psi(const F x[]) const; // get psi value at x
  F eval_psin(const F x[]) const; // get normalized psi at x

  bool magnetic_map(F rzp[3], F phi_end, int nsteps=100) const;
  bool magnetic_map_2pi_total_B(const ndarray<F>& totalB, F rzp[3], const int nsteps=3600) const;

  double theta(double r, double z) const;

  void derive_bfield0();
  void derive_curl_bfield0();

  ndarray<F> derive_delta_B(const ndarray<F>& As) const;
  ndarray<F> derive_total_B(const ndarray<F>& As) const;

public:
  std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> new_roi_mesh(
      std::vector<I> &node_map, /* vert id of original mesh --> vert id of new mesh */
      std::vector<I> &inverse_node_map /* vert id of new mesh --> vert id of original mesh */
  ) const; 

protected:
  xgc_units_t units;
  ndarray<F> psifield, bfield, bfield0; // b0 = B/|B|
  ndarray<F> curl_bfield0;
  ndarray<F> etemp_par, etemp_per, Te1d;

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
void simplicial_xgc_2d_mesh<I, F>::read_oneddiag(const std::string& filename, diy::mpi::communicator comm)
{
  try {
    etemp_par.read_file(filename, "/e_parallel_mean_en_avg");
    etemp_per.read_file(filename, "/e_perp_temperature_avg");
    Te1d = (etemp_par + etemp_per) * (2.0 / 3.0);
  } catch (...) {
    warn("unable to read e_parallel_mean_en_avg or e_perp_temperature_avg");
  }
}

template <typename I, typename F>
bool simplicial_xgc_2d_mesh<I, F>::eval_total_B(const ndarray<F>& totalB, const F rzp[3], F b[3]) const
{
  F mu[3];
  I tid = this->locator->locate(rzp, mu);
  if (tid < 0) {
    b[0] = b[1] = b[2] = F(0); 
    return false;
  } else {
    const F phin = mod2pi(rzp[2]) / (M_PI * 2);
    const int np = totalB.dim(2);
    const int i0 = (int(phin * np)) % np;
    const int i1 = (i0 + 1) % np;

    const F beta = phin * np - i0, alpha = F(1) - beta;

    F B[3][3];
    I tri[3];
    this->get_simplex(2, tid, tri);

    F rzp0[3], rzp1[3];
    for (int k = 0; k < 3; k ++) {
      rzp0[k] = rzp[k];
      rzp1[k] = rzp[k];
    }

    const int nsteps = 512; // TODO

    bool b0 = magnetic_map(rzp0, i0 * 2 * M_PI / np, beta * nsteps); 
    F B0[3] = {0}, B1[3] = {0};
    if (b0) {
      F mu[3];
      I tid = this->locate(rzp0, mu);
      if (tid < 0) 
        return false;
      else {
        I tri[3];
        this->get_simplex(2, tid, tri);
        for (int i = 0; i < 3; i ++) 
          for (int j = 0; j < 3; j ++)
            B0[j] += mu[i] * totalB(j, tri[i], i0);
      }
    } else 
      return false;
    bool b1 = magnetic_map(rzp1, i1 * 2 * M_PI / np, alpha * nsteps);
    if (b1) {
      F mu[3];
      I tid = this->locate(rzp1, mu);
      if (tid < 0) 
        return false;
      else {
        I tri[3];
        this->get_simplex(2, tid, tri);
        for (int i = 0; i < 3; i ++) 
          for (int j = 0; j < 3; j ++)
            B1[j] += mu[i] * totalB(j, tri[i], i1);
      }
    } else 
      return false;

    for (int i = 0; i < 3; i ++)
      b[i] = alpha * B0[i] + beta * B1[i];

    // FIXME: use magnetic map instead
#if 0
    for (int i = 0; i < 3; i ++) 
      for (int j = 0; j < 3; j ++)
        B[i][j] = alpha * totalB(j, tri[i], i0) + beta * totalB(j, tri[i], i1); //  bfield(j, tri[i]);
    lerp_s2v3(B, mu, b);
#endif
    
    // print3x3("B", B);
    // fprintf(stderr, "mu=%f, %f, %f, b=%f, %f, %f\n", 
    //     mu[0], mu[1], mu[2], b[0], b[1], b[2]);
    return true;
  }
}

template <typename I, typename F>
bool simplicial_xgc_2d_mesh<I, F>::magnetic_map_2pi_total_B(const ndarray<F>& totalB, F rzp[3], const int nsteps) const
{
  const F delta = 2 * M_PI / nsteps;

  for (int k = 0; k < nsteps; k ++) {
    if (!rk4<3, F>(rzp, [&](const F* rzp, F* v) {
          F B[3];
          if (eval_total_B(totalB, rzp, B)) {
            v[0] = rzp[0] * B[0] / B[2];
            v[1] = rzp[0] * B[1] / B[2];
            v[2] = 1;
            return true;
          } else 
            return false;
        }, delta))
      return false;
  }
  // rzp[2] = 0;
  return true;
}

template <typename I, typename F>
ndarray<F> simplicial_xgc_2d_mesh<I, F>::derive_total_B(const ndarray<F>& As) const
{
  ndarray<F> totalB;
  if (As.nd() == 1) {
    totalB = bfield + derive_delta_B(As);
    totalB.set_multicomponents();
  } else if (As.nd() == 2) {
    const auto deltaB = derive_delta_B(As);
    totalB.reshape(deltaB);
    for (int p = 0; p < As.dim(1); p ++)
      for (int i = 0; i < this->n(0); i ++)
        for (int j = 0; j < 3; j ++)
          totalB(j, i, p) = deltaB(j, i, p) + bfield(j, i);
    totalB.set_multicomponents(2);
  }
  return totalB;
}

template <typename I, typename F>
ndarray<F> simplicial_xgc_2d_mesh<I, F>::derive_delta_B(const ndarray<F>& As) const
{
  ndarray<F> deltaB;
  if (As.nd() == 1) { // wrong..
#if 0
    ndarray<F> gradAs = this->scalar_gradient(As);
    deltaB.reshape(3, this->n(0));
    deltaB.set_multicomponents();

    for (int i = 0; i < this->n(0); i ++) {
      deltaB(0, i) =  gradAs(1, i) * bfield0(2, i) - gradAs(2, i) * bfield0(1, i)  + As(i) * curl_bfield0(0, i);
      deltaB(1, i) = -gradAs(0, i) * bfield0(2, i) + gradAs(2, i) * bfield0(0, i)  + As(i) * curl_bfield0(1, i);
      deltaB(2, i) =  gradAs(0, i) * bfield0(1, i) - gradAs(1, i) * bfield0(0, i)  + As(i) * curl_bfield0(2, i);
    }
#endif
    assert(false);
  } else if (As.nd() == 2) {
    const int np = As.dim(1);
    const F dphi = 2 * M_PI / np;

    // ndarray<F> gradAs = this->vector_gradient(As.get_transpose()); // we only have dB/dR and dB/dZ

    deltaB.reshape(3, this->n(0), As.dim(1));
    deltaB.set_multicomponents(2);

    for (int c = 0; c < As.dim(1); c ++) {
      const int cnext = (c + 1) % np, 
                cprev = (c + np - 1) % np;

      ndarray<F> Asp;
      Asp.reshape(this->n(0));
      for (int i = 0; i < this->n(0); i ++)
        Asp[i] = As(i, c);
      ndarray<F> gradAs = this->scalar_gradient(As); // FIXME: this won't work as expected.  need to eval grad plane by plane

      for (int i = 0; i < this->n(0); i ++) {
        const F dAsdphi = (As(i, cnext) - As(i, cprev)) / (dphi * 2);

        // deltaB(0, i, c) =  gradAs(1, c, i) * bfield0(2, i) - gradAs(2, c, i) * bfield0(1, i)  + As(i, c) * curl_bfield0(0, i);
        // deltaB(1, i, c) = -gradAs(0, c, i) * bfield0(2, i) + gradAs(2, c, i) * bfield0(0, i)  + As(i, c) * curl_bfield0(1, i);
        // deltaB(2, i, c) =  gradAs(0, c, i) * bfield0(1, i) - gradAs(1, c, i) * bfield0(0, i)  + As(i, c) * curl_bfield0(2, i);
        deltaB(0, i, c) =  gradAs(1, i) * bfield0(2, i) - dAsdphi      * bfield0(1, i)  + As(i, c) * curl_bfield0(0, i);
        deltaB(1, i, c) = -gradAs(0, i) * bfield0(2, i) + dAsdphi      * bfield0(0, i)  + As(i, c) * curl_bfield0(1, i);
        deltaB(2, i, c) =  gradAs(0, i) * bfield0(1, i) - gradAs(1, i) * bfield0(0, i)  + As(i, c) * curl_bfield0(2, i);
      }
    }
  }

  return deltaB;
}

template <typename I, typename F>
void simplicial_xgc_2d_mesh<I, F>::derive_curl_bfield0()
{
  if (bfield0.empty()) derive_bfield0();
  
  ndarray<F> gradRZ_bfield0 = this->vector_gradient(bfield0); // dbr/dphi, dbz/dphi, dbphi/dphi are zero
  curl_bfield0.reshape(3, bfield0.dim(1));

  for (int i = 0; i < bfield0.dim(1); i ++) {
    const F R = this->vertex_coords(0, i);
    curl_bfield0(0, i) = -gradRZ_bfield0(1, 2, i); // R
    curl_bfield0(1, i) = bfield0(2, i) / R + gradRZ_bfield0(0, 2, i); // Z
    curl_bfield0(2, i) = gradRZ_bfield0(1, 0, i) - gradRZ_bfield0(0, 1, i); // Phi
  }
}

template <typename I, typename F>
void simplicial_xgc_2d_mesh<I, F>::derive_bfield0()
{
  bfield0.reshape(bfield);
  bfield0.set_multicomponents();
  for (int i = 0; i < bfield.dim(1); i ++) {
    const F br = bfield(0, i), 
            bz = bfield(1, i),
            bphi = bfield(2, i);
    const F norm = std::sqrt(br*br + bz*bz + bphi*bphi);
    bfield0(0, i) = br / norm;
    bfield0(1, i) = bz / norm;
    bfield0(2, i) = bphi / norm;
  }
}

template <typename I, typename F>
void simplicial_xgc_2d_mesh<I, F>::read_bfield(const std::string& filename, diy::mpi::communicator comm)
{
  auto ext = file_extension(filename);
  if (ext == FILE_EXT_BP) bfield.read_bp(filename, "/node_data[0]/values", -1, comm);
  else if (ext == FILE_EXT_HDF5) bfield.read_h5(filename, "/node_data[0]/values");
  else fatal(FTK_ERR_FILE_UNRECOGNIZED_EXTENSION);

  bfield.set_multicomponents();
  // derive_bfield0();
  derive_curl_bfield0();
}

template <typename I, typename F>
bool simplicial_xgc_2d_mesh<I, F>::magnetic_map(F rzp[3], F phi_end, int nsteps) const
{
  if (nsteps == 0) return true;
  if (bfield.empty()) 
    fatal("missing xgc bfield.");

  const F delta = (phi_end - rzp[2]) / nsteps;

  for (int k = 0; k < nsteps; k ++) {
    if (!rk4<3, F>(rzp, [&](const F* rzp, F* v) {
          F B[3];
          if (eval_b(rzp, B)) {
            v[0] = rzp[0] * B[0] / B[2];
            v[1] = rzp[0] * B[1] / B[2];
            v[2] = 1;
            return true;
          } else 
            return false;
        }, delta))
      return false;

    // fprintf(stderr, "integrating %f, %f, %f\n", rzp[0], rzp[1], rzp[2]);
  }
  return true;
}

template <typename I, typename F>
void simplicial_xgc_2d_mesh<I, F>::initialize_point_locator()
{
  this->locator.reset( new point_locator_2d_quad<I, F>(*this) );
}

template <typename I, typename F>
F simplicial_xgc_2d_mesh<I, F>::eval_psi(const F x[]) const
{
  F val;
  if (simplicial_unstructured_2d_mesh<I, F>::eval(psifield, x, &val))
    return val;
  else return NAN;
}

template <typename I, typename F>
F simplicial_xgc_2d_mesh<I, F>::eval_psin(const F x[]) const
{
  return eval_psi(x) / units.psi_x;
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
std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> simplicial_xgc_2d_mesh<I, F>::from_xgc_mesh_file(
    const std::string& filename,
    diy::mpi::communicator comm)
{
  auto ext = file_extension(filename);
  if (ext == FILE_EXT_BP) return from_xgc_mesh_bp(filename, comm);
  else if (ext == FILE_EXT_HDF5) return from_xgc_mesh_h5(filename);
  else {
    fatal(FTK_ERR_FILE_UNRECOGNIZED_EXTENSION);
    return NULL;
  }
}

template <typename I, typename F>
std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> simplicial_xgc_2d_mesh<I, F>::from_xgc_mesh_bp(
    const std::string& filename,
    diy::mpi::communicator comm)
{
  ndarray<I> triangles = ndarray<I>::from_bp(filename, "/cell_set[0]/node_connect_list", -1, comm);
  ndarray<F> coords = ndarray<F>::from_bp(filename, "/coordinates/values", -1, comm);
  ndarray<I> nextnodes = ndarray<I>::from_bp(filename, "nextnode", -1, comm);
  ndarray<F> psifield = ndarray<F>::from_bp(filename, "psi", -1, comm);

  // triangles.transpose();
  // coords.transpose();

  // std::cerr << triangles.shape() << std::endl;
  // std::cerr << coords.shape() << std::endl;
  // std::cerr << nextnodes.shape() << std::endl;
  // std::cerr << psifield.shape() << std::endl;

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
std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> simplicial_xgc_2d_mesh<I, F>::new_roi_mesh(
    std::vector<I> &node_map, 
    std::vector<I> &inverse_node_map) const
{
  // if (roi_nodes.empty()) 
  //   initialize_roi();
  const auto nn = roi_nodes.size();

  std::map<I, I> map; // oridinal id to new id
  ndarray<F> new_coords({2, nn});
  for (auto i = 0; i < nn; i ++)
    map[roi_nodes[i]] = i;

  for (auto i = 0; i < nn; i ++) {
    I node = roi_nodes[i];
    for (int j = 0; j < 2; j ++)
      new_coords(j, i) = this->vertex_coords(j, node);
  }

  // updating node map
  node_map.resize(this->n(0), -1);
  inverse_node_map.resize(nn);
  for (auto i = 0; i < this->n(0); i ++) {
    if (map.find(i) != map.end()) {
      node_map[i] = map[i];
      inverse_node_map[map[i]] = i;
    }
  }
 
  // new triangles
  std::vector<I> new_triangles;
  for (auto i = 0; i < this->n(2); i ++) {
    I tri[3];
    this->get_triangle(i, tri);

    bool b = true;
    for (int j = 0; j < 3; j ++)
      if (node_map[tri[j]] < 0)
        b = false;
    if (!b) continue;  // not a triangle in the new mesh

    for (int j = 0; j < 3; j ++)
      new_triangles.push_back( map[tri[j]] );
  }

  ndarray<I> new_triangles_array;
  new_triangles_array.copy_vector(new_triangles);
  new_triangles_array.reshape(3, new_triangles.size() / 3);
  
  // std::shared_ptr<simplicial_unstructured_2d_mesh<I, F>> m2( new simplicial_unstructured_2d_mesh<I, F>( new_coords, new_triangles_array ) );
  // return m2;

  std::shared_ptr<simplicial_xgc_2d_mesh<I, F>> mx2( new simplicial_xgc_2d_mesh<I, F>( new_coords, new_triangles_array ) );
  mx2->units = units;
 
  // updating psifield and bfield
  if (!psifield.empty()) {
    mx2->psifield.reshape(nn);
    for (int i = 0; i < nn; i ++) {
      const auto k = roi_nodes[i];
      mx2->psifield[i] = psifield[k];
    }
  }

  if (!bfield.empty()) {
    mx2->bfield.reshape(3, nn);
    for (int i = 0; i < nn; i ++) {
      const auto k = roi_nodes[i];
      mx2->psifield[i] = psifield[k];
      for (int j = 0; j < 3; j ++) 
        mx2->bfield(j, i) = bfield(j, k);
    }
  }
  return mx2;
}

template <typename I, typename F>
double simplicial_xgc_2d_mesh<I, F>::theta(double r, double z) const 
{
  return atan2(z - units.eq_axis_z, r - units.eq_axis_r);
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
