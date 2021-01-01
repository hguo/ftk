#ifndef _FTK_XGC_BLOB_CCL_TRACKER_HH
#define _FTK_XGC_BLOB_CCL_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/filters/xgc_tracker.hh>
#include <ftk/tracking_graph/tracking_graph.hh>

namespace ftk {

struct xgc_blob_threshold_tracker : public xgc_tracker {
  xgc_blob_threshold_tracker(diy::mpi::communicator comm, 
      std::shared_ptr<simplicial_xgc_3d_mesh<>> m3) : xgc_tracker(comm, m3) {}
  virtual ~xgc_blob_threshold_tracker() {}
  
  int cpdims() const { return 0; }
  
  void initialize() {}
  void update() {}
  void finalize() {}

  void update_timestep();

  void push_field_data_snapshot(const ndarray<double> &scalar) {}
};

/////

inline void xgc_blob_threshold_tracker::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);

  // m3->element_for(0, current_timestep, func, xl, nthreads, enable_set_affinity);
}

}

#endif
