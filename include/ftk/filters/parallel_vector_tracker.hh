#ifndef _FTK_PV_TRACKER_HH
#define _FTK_PV_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/filters/feature_point.hh>
#include <ftk/filters/feature_surface.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/external/diy/serialization.hpp>
#include <ftk/external/diy-ext/gather.hh>

namespace ftk {

struct parallel_vector_tracker : public virtual tracker {
  parallel_vector_tracker(diy::mpi::communicator comm) : tracker(comm) {}

  void reset() { field_data_snapshots.clear(); }

public:
  bool advance_timestep();

public:
  bool pop_field_data_snapshot();
  void push_field_data_snapshot(
      const ndarray<double> &v, 
      const ndarray<double> &w, 
      const ndarray<double> &Jv);
  
protected:
  struct field_data_snapshot_t {
    ndarray<double> v,  w;
    ndarray<double> Jv; // jacobian of v
  };
  std::deque<field_data_snapshot_t> field_data_snapshots;
};

////////
inline bool parallel_vector_tracker::advance_timestep()
{
  update_timestep();
  pop_field_data_snapshot();

  current_timestep ++;
  return field_data_snapshots.size() > 0;
}

inline bool parallel_vector_tracker::pop_field_data_snapshot()
{
  if (field_data_snapshots.size() > 0) {
    field_data_snapshots.pop_front();
    return true;
  } else return false;
}

inline void parallel_vector_tracker::push_field_data_snapshot(
    const ndarray<double>& v,
    const ndarray<double>& w, 
    const ndarray<double>& Jv)
{
  field_data_snapshot_t snapshot;
  snapshot.v = v;
  snapshot.w = w;
  snapshot.Jv = Jv;

  field_data_snapshots.emplace_back(snapshot);
}


}

#endif
