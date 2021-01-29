#ifndef _FTK_XGC_TRACKER_HH
#define _FTK_XGC_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/filters/feature_point.hh>
#include <ftk/filters/feature_surface.hh>
#include <ftk/filters/feature_volume.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/utils/gather.hh>
#include <ftk/mesh/simplicial_xgc_2d_mesh.hh>
#include <ftk/mesh/simplicial_xgc_3d_mesh.hh>
#include <ftk/mesh/simplicial_xgc_3dff_mesh.hh>
#include <iomanip>

namespace ftk {

struct xgc_tracker : public tracker {
  xgc_tracker(diy::mpi::communicator comm,
      std::shared_ptr<simplicial_xgc_3d_mesh<>> mx);

  void reset() { field_data_snapshots.clear(); }

public:
  std::shared_ptr<simplicial_unstructured_2d_mesh<>> get_m2() { return m2; }
  void initialize_ff_mesh(const std::string& filename) { mf3.reset(new simplicial_xgc_3dff_mesh<>(m2, m3->get_nphi(), m3->get_iphi(), m3->get_vphi()) ); mf3->initialize_ff_mesh(filename); }

public:
  virtual void push_field_data_snapshot(const ndarray<double> &scalar);
  virtual void push_field_data_snapshot(
      const ndarray<double> &scalar, 
      const ndarray<double> &vector,
      const ndarray<double> &jacobian);
  
  bool pop_field_data_snapshot();

  bool advance_timestep();

protected:
  struct field_data_snapshot_t {
    ndarray<double> scalar, vector, jacobian;
  };
  std::deque<field_data_snapshot_t> field_data_snapshots;

  std::shared_ptr<simplicial_xgc_2d_mesh<>> m2; 
  std::shared_ptr<simplicial_xgc_3d_mesh<>> m3; 
  std::shared_ptr<simplicial_unstructured_extruded_3d_mesh<>> m4;
  
  std::shared_ptr<simplicial_xgc_3dff_mesh<>> mf3; // field following 3d mesh
};

//////
xgc_tracker::xgc_tracker(
    diy::mpi::communicator comm, 
    std::shared_ptr<simplicial_xgc_3d_mesh<>> mx) :
  tracker(comm),
  m2(mx->get_m2()),
  m3(mx),
  m4(new ftk::simplicial_unstructured_extruded_3d_mesh<>(*m3))
{

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
  ndarray<double> F, G, J;

  F.reshape(scalar);
  G.reshape(2, scalar.dim(0), scalar.dim(1));
  J.reshape(2, 2, scalar.dim(0), scalar.dim(1));
  
  for (size_t i = 0; i < m3->get_nphi(); i ++) {
    ftk::ndarray<double> f, grad, j;
    auto slice = scalar.slice_time(i);
    m2->smooth_scalar_gradient_jacobian(slice, f, grad, j);
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

}

#endif
