#ifndef _FTK_XGC_BLOB_FILAMENT_TRACKER_HH
#define _FTK_XGC_BLOB_FILAMENT_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/filters/feature_point.hh>
#include <ftk/filters/feature_surface.hh>
#include <ftk/filters/feature_volume.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/external/diy/serialization.hpp>
#include <ftk/external/diy-ext/gather.hh>
#include <ftk/io/tdgl.hh>
#include <iomanip>

namespace ftk {
  
struct xgc_blob_filament_tracker : public virtual tracker {
  xgc_blob_filament_tracker(diy::mpi::communicator comm, 
      const simplicial_unstructured_2d_mesh<>& m2, 
      int nphi_, int iphi_);

  int cpdims() const { return 3; }
 
  void set_nphi_iphi(int nphi, int iphi);

  void initialize() {}
  void update();
  void reset() {
    field_data_snapshots.clear();
  }
  void finalize() {}
  
public:
  void update_timestep() {}
  bool advance_timestep();

public:
  bool pop_field_data_snapshot();
  void push_field_data_snapshot(
      const ndarray<double> &scalar, 
      const ndarray<double> &vector,
      const ndarray<double> &jacobian);
  
protected:
  struct field_data_snapshot_t {
    ndarray<double> scalar, vector, jacobian;
  };
  std::deque<field_data_snapshot_t> field_data_snapshots;

  const int nphi, iphi;
  const simplicial_unstructured_2d_mesh<>& m2;
  const simplicial_unstructured_3d_mesh<> m3;
  const simplicial_unstructured_extruded_3d_mesh<> m4;
};

/////
  
xgc_blob_filament_tracker::xgc_blob_filament_tracker(
    diy::mpi::communicator comm, 
    const simplicial_unstructured_2d_mesh<>& m2_, 
    int nphi_, int iphi_) :
  tracker(comm),
  m2(m2_), nphi(nphi_), iphi(iphi_),
  m3(ftk::simplicial_unstructured_3d_mesh<>::from_xgc_mesh(m2_, nphi_, iphi_)),
  m4(m3)
{

}

inline bool xgc_blob_filament_tracker::advance_timestep()
{
  update_timestep();
  pop_field_data_snapshot();

  current_timestep ++;
  return field_data_snapshots.size() > 0;
}

inline bool xgc_blob_filament_tracker::pop_field_data_snapshot()
{
  if (field_data_snapshots.size() > 0) {
    field_data_snapshots.pop_front();
    return true;
  } else return false;
}
  
inline void xgc_blob_filament_tracker::push_field_data_snapshot(
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

inline void xgc_blob_filament_tracker::update()
{
  // WIP
}

}

#endif
