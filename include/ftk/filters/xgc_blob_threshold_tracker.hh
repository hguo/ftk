#ifndef _FTK_XGC_BLOB_CCL_TRACKER_HH
#define _FTK_XGC_BLOB_CCL_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/tracking_graph/tracking_graph.hh>

namespace ftk {

struct xgc_blob_threshold_tracker : public connected_component_tracker<int, int> {
  xgc_blob_threshold_tracker(diy::mpi::communicator comm, 
      std::shared_ptr<simplicial_xgc_3d_mesh> m3);
  virtual ~xgc_blob_threshold_tracker() {}
  
  void update_timestep();
  bool advance_timestep();

  void push_field_data_snapshot(const ndarray<double> &scalar);
  bool pop_field_data_snapshot();

protected:
  struct field_data_snapshot_t {
    ndarray<double> scalar;
  };

protected:
  std::shared_ptr<simplicial_xgc_2d_mesh> m2;
  std::shared_ptr<simplicial_xgc_3d_mesh> m3;
};

/////

inline xgc_blob_threshold_tracker::xgc_blob_threshold_tracker(
    diy::mpi::communicator comm, 
    std::shared_ptr<simplicial_xgc_3d_mesh> mx) :
  m2(mx), m3(mx->get_m2())
{
}

inline void xgc_blob_threshold_tracker::push_field_data_snapshot(
    const ndarray<double> &scalar)
{
  field_data_snapshot_t snapshot;
  snapshot.scalar = scalar;

  field_data_snapshots.emplace_back(snapshot);
}

inline bool xgc_blob_filament_tracker::advance_timestep()
{
  update_timestep();
  pop_field_data_snapshot();

  current_timestep ++;
  return field_data_snapshots.size() > 0;
}

inline void xgc_blob_threshold_tracker::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);

  m4->element_for_ordinal(0, current_timestep, func, xl, nthreads, enable_set_affinity);
}

}

#endif
