#ifndef _FTK_XGC_TRACKER_HH
#define _FTK_XGC_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/features/feature_point.hh>
#include <ftk/features/feature_surface.hh>
#include <ftk/features/feature_volume.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/utils/gather.hh>
#include <ftk/ndarray/field_data_snapshot.hh>
#include <ftk/mesh/simplicial_xgc_2d_mesh.hh>
#include <ftk/mesh/simplicial_xgc_3d_mesh.hh>
#include <ftk/mesh/simplicial_xgc_3dff_mesh.hh>
#include <ftk/mesh/simplicial_unstructured_extruded_3d_mesh.hh>
#include <iomanip>

namespace ftk {

struct xgc_tracker : public tracker {
  xgc_tracker(diy::mpi::communicator comm,
      std::shared_ptr<simplicial_xgc_3d_mesh<>> mx);

  void reset() { field_data_snapshots.clear(); }

  void set_use_roi(bool b) { use_roi = b; }

public:
  std::shared_ptr<simplicial_unstructured_2d_mesh<>> get_m2() { return m2; }
  void initialize_ff_mesh(const std::string& filename) { mf3.reset(new simplicial_xgc_3dff_mesh<>(m2, m3->get_nphi(), m3->get_iphi(), m3->get_vphi()) ); mf3->initialize_ff_mesh(filename); }

public:
  virtual void push_field_data_snapshot(std::shared_ptr<ndarray_group>);
  virtual void push_field_data_snapshot(const ndarray<double> &scalar);
  virtual void push_field_data_snapshot(
      const ndarray<double> &scalar, 
      const ndarray<double> &vector,
      const ndarray<double> &jacobian);
  
  bool pop_field_data_snapshot();

  bool advance_timestep();

public:
  double derive_threshold(const ndarray<double>& scalar, double sigma_threshold = 2.5) const;

protected:
  struct field_data_snapshot_t {
    ndarray<double> scalar, vector, jacobian;
    ndarray<double> Er;
  };
  std::deque<field_data_snapshot_t> field_data_snapshots;

  std::shared_ptr<simplicial_xgc_2d_mesh<>> m2; 
  std::shared_ptr<simplicial_xgc_3d_mesh<>> m3, m30;
  std::shared_ptr<simplicial_unstructured_extruded_3d_mesh<>> m4;

  // roi meshes
  bool use_roi = false;
  std::shared_ptr<simplicial_xgc_2d_mesh<>> mr2; 
  std::shared_ptr<simplicial_xgc_3d_mesh<>> mr3; 
  std::shared_ptr<simplicial_unstructured_extruded_3d_mesh<>> mr4;
  std::vector<int> roi_node_map, roi_inverse_node_map;

  std::shared_ptr<simplicial_xgc_3dff_mesh<>> mf3; // field following 3d mesh
};

//////
xgc_tracker::xgc_tracker(
    diy::mpi::communicator comm, 
    std::shared_ptr<simplicial_xgc_3d_mesh<>> mx) :
  tracker(comm),
  m2(mx->get_m2()),
  m3(mx),
  m4(new simplicial_unstructured_extruded_3d_mesh<>(*m3))
{
  // initialize roi meshes
  m2->initialize_roi();
  mr2 = m2->new_roi_mesh(roi_node_map, roi_inverse_node_map);
  mr3.reset(new simplicial_xgc_3d_mesh<>(mr2, 
        m3->get_nphi(), m3->get_iphi(), m3->get_vphi()));
  m30.reset(new simplicial_xgc_3d_mesh<>(m2, 
        m3->get_nphi(), m3->get_iphi()));
  mr4.reset(new simplicial_unstructured_extruded_3d_mesh<>(*mr3));
}

inline void xgc_tracker::push_field_data_snapshot(std::shared_ptr<ndarray_group> g) 
{
  field_data_snapshot_t s;
  
  auto density = g->get<double>("density"); // .get_transpose();
  if (density.dim(0) < density.dim(1))
    density.transpose();
  m30->smooth_scalar_gradient_jacobian(density, s.scalar, s.vector, s.jacobian);
  
  if (g->has("Er")) {
    // auto Er = 
    // s.Er = m30->smooth_scalar(Er);
    s.Er = g->get<double>("Er"); // .get_transpose();
    if (s.Er.dim(0) < s.Er.dim(1))
      s.Er.transpose();
  }
  
  field_data_snapshots.emplace_back(s);
}

inline void xgc_tracker::push_field_data_snapshot(
      const ndarray<double> &scalar, 
      const ndarray<double> &vector,
      const ndarray<double> &jacobian)
{
  field_data_snapshot_t snapshot;
  snapshot.scalar = scalar;
  snapshot.vector = vector;
  snapshot.jacobian = jacobian;

  field_data_snapshots.emplace_back(snapshot);
}

inline void xgc_tracker::push_field_data_snapshot(
      const ndarray<double> &scalar)
{
#if 0
  auto filename = series_filename("xgc-smoothed-%05d.vtu", current_timestep);
  m30->scalar_to_vtu_slices_file(filename, "scalar", scalar);
  return;
  // m30->scalar_to_vtu_slices_file(filename, "vector", G);
#endif

  ndarray<double> F, G, J;

  F.reshape(scalar);
  G.reshape(2, scalar.dim(0), scalar.dim(1));
  G.set_multicomponents(1);
  J.reshape(2, 2, scalar.dim(0), scalar.dim(1));
  J.set_multicomponents(2);
  
  for (size_t i = 0; i < m3->get_nphi(); i ++) {
    ndarray<double> f, grad, j;
    auto slice = scalar.slice_time(i);
    m2->smooth_scalar_gradient_jacobian(slice, f, grad, j);
    // f = m2->smooth_scalar(slice);
    // m2->scalar_gradient_jacobian(slice, grad, j);
    // m2->scalar_gradient_jacobian(f, grad, j);
    // if (i == 0)
    //   m2->array_to_vtu("nonsmooth-grad.vtu", "grad", grad);
    // f = slice;
    for (size_t k = 0; k < m2->n(0); k ++) {
      F(k, i) = f(k);
      G(0, k, i) = grad(0, k);
      G(1, k, i) = grad(1, k);
      J(0, 0, k, i) = j(0, 0, k);
      J(1, 0, k, i) = j(1, 0, k);
      J(1, 1, k, i) = j(1, 1, k);
      J(0, 1, k, i) = j(0, 1, k);
    }
  }

  push_field_data_snapshot(F, G, J);
}
  
inline bool xgc_tracker::pop_field_data_snapshot()
{
  if (field_data_snapshots.size() > 0) {
    field_data_snapshots.pop_front();
    return true;
  } else return false;
}

inline bool xgc_tracker::advance_timestep()
{
  update_timestep();
  pop_field_data_snapshot();

  current_timestep ++;
  return field_data_snapshots.size() > 0;

#if 0
  // fprintf(stderr, "advancing timestep!!! t=%d, #snapshot=%zu\n", 
  //     current_timestep, field_data_snapshots.size());

  if (field_data_snapshots.size() >= 2) {
    update_timestep();
    pop_field_data_snapshot();
    current_timestep ++;
  }

  return field_data_snapshots.size() > 0;
#endif
}

inline double xgc_tracker::derive_threshold(const ndarray<double>& scalar, double sigma_threshold) const
{
  const std::vector<int> roi_nodes = m2->get_roi_nodes();
  double sum = 0.0, sigma = 0.0;
  const int nphi = m3->get_nphi(), m2n0 = m2->n(0), nroi0 = roi_nodes.size();

  for (int p = 0; p < nphi; p ++) {
    for (int i = 0; i < nroi0; i ++) {
      const int idx = roi_nodes[i] + m2n0 * p;
      sum += scalar[idx];
    }
  }

  const double ave = sum / nroi0;
  for (int p = 0; p < nphi; p ++) {
    for (int i = 0; i < nroi0; i ++) {
      const int idx = roi_nodes[i] + m2n0 * p;
      sigma += std::pow(scalar[idx] - ave, 2.0);
    }
  }

  sigma = std::sqrt(sigma / (nphi * nroi0)) * sigma_threshold;
  return sigma;
}

}

#endif
