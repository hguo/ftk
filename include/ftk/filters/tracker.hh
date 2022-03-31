#ifndef _FTK_TRACKER_HH
#define _FTK_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/ndarray/field_data_snapshot.hh>
#include <ftk/ndarray/ndarray_group.hh>
#include <ftk/filters/filter.hh>
#include <ftk/external/diy/master.hpp>

namespace ftk {

enum {
  TRACKER_CRITICAL_POINT = 1,
  TRACKER_CRITICAL_LINE = 2,
  TRACKER_SUJUDI_HAIMES = 2,
  TRACKER_LEVY_DEGANI_SEGINER = 2,
  TRACKER_RIDGE_VALLEY = 2,
  TRACKER_TDGL_VORTEX = 3,
  TRACKER_CONTOUR = 4,
  TRACKER_CONNECTED_COMPONENTS = 5,
  TRACKER_THRESHOLD = 6,
  TRACKER_PARTICLE = 7,
  TRACKER_XGC_BLOB_FILAMENT = 105,
  TRACKER_XGC_BLOB_THRESHOLD = 106,
  TRACKER_MPAS_O_CRITICAL_POINT = 201
};

struct tracker : public filter
{
  tracker(diy::mpi::communicator comm) : filter(comm) {} // , master(comm) {}
  virtual ~tracker() {}
  
  // virtual int cpdims() const = 0; // featutre dimension
  
  void set_start_timestep(int t) { start_timestep = t;}
  void set_end_timestep(int t) { end_timestep = t; }
  
  virtual void set_current_timestep(int t) {current_timestep = t;}
  int get_current_timestep() const {return current_timestep;}
 
  void set_input_array_partial(bool b) {is_input_array_partial = b;}
  void set_use_default_domain_partition(bool b) {use_default_domain_partition = true;}

  static int str2tracker(const std::string&);

public:
  virtual void initialize() = 0;
  virtual void finalize() = 0;
  
  virtual bool advance_timestep() = 0;
  virtual void update_timestep() = 0;

public:
  virtual void push_field_data_snapshot(std::shared_ptr<ndarray_group> g) {snapshots.push_back(g);}
  virtual void push_field_data_snapshot(const std::string key, const ndarray<double>& arr) {
    std::shared_ptr<ndarray_group> g(new ndarray_group);
    g->set(key, arr);
    push_field_data_snapshot(g);
  }
  virtual bool pop_field_data_snapshot();

protected:
  std::deque<std::shared_ptr<ndarray_group>> snapshots;

protected:
  // diy::Master master;
  
protected:
  int start_timestep = 0, 
      end_timestep = std::numeric_limits<int>::max();

  int current_timestep = 0;
 
protected:
  bool is_input_array_partial = false;
  bool use_default_domain_partition = true;

protected: // benchmark
  double accumulated_kernel_time = 0.0;
};

////////
inline int tracker::str2tracker(const std::string& s) 
{
  if (s == "cp" || s == "critical_point") 
    return TRACKER_CRITICAL_POINT;
  else if (s == "iso" || s == "isovolume" || s == "isosurface" || s == "isosurfaces")
    return TRACKER_CONTOUR;
  else if (s == "tdgl" || s == "tdgl_vortex" || s == "tdgl-vortex" || s == "tdgl_vortices" || s == "tdgl-vortices")
    return TRACKER_TDGL_VORTEX;
  else if (s == "cl" || s == "critical_line" || s == "critical_lines")
    return TRACKER_CRITICAL_LINE;
  else if (s == "mpas-o-cp")
    return TRACKER_MPAS_O_CRITICAL_POINT;
  else if (s == "sujudi_haimes")
    return TRACKER_SUJUDI_HAIMES;
  else if (s == "ridge_valley")
    return TRACKER_RIDGE_VALLEY;
  else if (s == "levy_degani_seginer")
    return TRACKER_LEVY_DEGANI_SEGINER;
  else if (s == "particle" || s == "pt")
    return TRACKER_PARTICLE;
  else if (s == "cc" || s == "connected_component" || s == "connected_components")
    return TRACKER_CONNECTED_COMPONENTS;
  else if (s == "xgc_blob_filament" || s == "xgc-blob-filament")
    return TRACKER_XGC_BLOB_FILAMENT;
  else if (s == "xgc_blob_threshold" || s == "xgc-blob-threshold")
    return TRACKER_XGC_BLOB_THRESHOLD;
  else return 0;
}

inline bool tracker::pop_field_data_snapshot()
{
  if (snapshots.size() > 0) {
    snapshots.pop_front();
    return true;
  } else return false;
}

}

#endif
