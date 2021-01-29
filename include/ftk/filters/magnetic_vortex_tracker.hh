#ifndef _FTK_MAGNETIC_VORTEX_HH
#define _FTK_MAGNETIC_VORTEX_HH

#include <ftk/config.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/filters/feature_point.hh>
#include <ftk/filters/feature_surface.hh>
#include <ftk/filters/feature_volume.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/utils/gather.hh>
#include <ftk/io/magnetic.hh>
#include <iomanip>

namespace ftk {

struct magnetic_vortex_tracker : public virtual tracker {
  magnetic_vortex_tracker(diy::mpi::communicator comm) : tracker(comm) {}

  int cpdims() const { return 3; }

  virtual void update() {}; 
  void reset() {
    field_data_snapshots.clear();
  }

public:
  bool advance_timestep();

public:
  bool pop_field_data_snapshot();
  void push_field_data_snapshot(const ndarray<float> &spin);

protected:
  struct field_data_snapshot_t {
    ndarray<float> spin
      std::deque<field_data_snapshot_t> field_data_snapshots;;
  };
  std::deque<field_data_snapshot_t> field_data_snapshots;
};

inline bool magnetic_vortex_tracker::advance_timestep()
{
  update_timestep();
  pop_field_data_snapshot();

  current_timestep ++;
  return field_data_snapshots.size() > 0;
}

inline bool magnetic_vortex_tracker::pop_field_data_snapshot()
{
  if (field_data_snapshots.size() > 0) {
    field_data_snapshots.pop_front();
    return true;
  } else return false;
}
  
inline void magnetic_vortex_tracker::push_field_data_snapshot(const ndarray<float> &spin)
{
  field_data_snapshot_t snapshot;
  snapshot.spin = spin;
  field_data_snapshots.emplace_back(snapshot);
}

}

#endif
