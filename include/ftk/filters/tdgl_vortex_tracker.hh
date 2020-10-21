#ifndef _FTK_TDGL_TRACKER_HH
#define _FTK_TDGL_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/filters/filter.hh>
#include <ftk/filters/feature_point.hh>
#include <ftk/filters/feature_surface.hh>
#include <ftk/filters/feature_volume.hh>
#include <ftk/filters/feature_volume_set.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/external/diy/serialization.hpp>
#include <ftk/external/diy-ext/gather.hh>
#include <ftk/io/tdgl.hh>
#include <iomanip>

namespace ftk {
  
struct tdgl_vortex_tracker : public filter {
  int cpdims() const { return 3; }
  
  virtual void update() {}; 
  void reset() {
    field_data_snapshots.clear();
    // traced_contours.clear();
  }
  
public:
  virtual void initialize() = 0;
  virtual void finalize() = 0;

  bool advance_timestep();
  virtual void update_timestep() = 0;

public:
  bool pop_field_data_snapshot();
  virtual void push_field_data_snapshot(
      const tdgl_metadata_t &meta,
      const ndarray<float> &rho_phi,
      const ndarray<float> &re_im
  );

  virtual void set_current_timestep(int t) {current_timestep = t;}
  int get_current_timestep() const {return current_timestep;}

  void set_end_timestep(int t) {end_timestep = t;}
  
protected:
  struct field_data_snapshot_t {
    tdgl_metadata_t meta;
    ndarray<float> re_im, rho_phi;
  };
  std::deque<field_data_snapshot_t> field_data_snapshots;
  
  int current_timestep = 0;
  int start_timestep = 0, 
      end_timestep = std::numeric_limits<int>::max();
};

/////
inline bool tdgl_vortex_tracker::advance_timestep()
{
  update_timestep();
  pop_field_data_snapshot();

  current_timestep ++;
  return field_data_snapshots.size() > 0;
}

}

#endif
