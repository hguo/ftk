#ifndef _FTK_TRACKER_HH
#define _FTK_TRACKER_HH

#include <ftk/config.hh>
#include <ndarray/ndarray_group.hh>
#include <ftk/filters/filter.hh>
#include <ftk/external/diy/master.hpp>

namespace ftk {

/**
 * @brief Tracker type identifiers
 *
 * These constants identify different types of feature trackers available in FTK.
 * Use with tracker::str2tracker() to convert string names to tracker types.
 */
enum {
  TRACKER_CRITICAL_POINT = 1,        ///< Critical point tracker for scalar/vector fields
  TRACKER_CRITICAL_LINE = 2,         ///< Critical line tracker for 3D vector fields
  TRACKER_SUJUDI_HAIMES = 2,         ///< Sujudi-Haimes critical line detection
  TRACKER_LEVY_DEGANI_SEGINER = 2,   ///< Levy-Degani-Seginer critical line detection
  TRACKER_RIDGE_VALLEY = 2,          ///< Ridge and valley detection
  TRACKER_TDGL_VORTEX = 3,           ///< TDGL vortex tracker
  TRACKER_CONTOUR = 4,               ///< Isosurface/contour tracker
  TRACKER_CONNECTED_COMPONENTS = 5,  ///< Connected component tracker
  TRACKER_THRESHOLD = 6,             ///< Threshold-based tracker
  TRACKER_PARTICLE = 7,              ///< Particle advection tracker
  TRACKER_XGC_BLOB_FILAMENT = 105,   ///< XGC blob filament tracker
  TRACKER_XGC_BLOB_THRESHOLD = 106,  ///< XGC blob threshold tracker
  TRACKER_MPAS_O_CRITICAL_POINT = 201, ///< MPAS-Ocean critical point tracker
  TRACKER_MPAS_O_PARTICLES = 202     ///< MPAS-Ocean particle tracker
};

/**
 * @brief Base class for all feature trackers
 *
 * The tracker class provides the foundation for tracking features across time
 * in scientific datasets. It manages time-series data through snapshots and
 * provides a streaming interface for processing temporal sequences.
 *
 * Trackers operate in a three-phase workflow:
 * 1. initialize() - Set up the tracker with domain and configuration
 * 2. Push snapshots and update_timestep() - Process each timestep incrementally
 * 3. finalize() - Complete tracking and produce final results
 *
 * Derived classes implement specific tracking algorithms for different feature
 * types (critical points, contours, particles, etc.).
 */
struct tracker : public filter
{
  tracker(diy::mpi::communicator comm) : filter(comm) {} // , master(comm) {}
  virtual ~tracker() {}

  // virtual int cpdims() const = 0; // featutre dimension

  /**
   * @brief Set the starting timestep for tracking
   * @param t Starting timestep index
   */
  void set_start_timestep(int t) { start_timestep = t;}

  /**
   * @brief Set the ending timestep for tracking
   * @param t Ending timestep index
   */
  void set_end_timestep(int t) { end_timestep = t; }

  /**
   * @brief Set the total number of timesteps
   * @param n Total number of timesteps (0 for unlimited)
   */
  void set_ntimesteps(int n) { ntimesteps = n; }

  /**
   * @brief Set the current timestep
   * @param t Current timestep index
   */
  virtual void set_current_timestep(int t) {current_timestep = t;}

  /**
   * @brief Get the current timestep
   * @return Current timestep index
   */
  int get_current_timestep() const {return current_timestep;}

  /**
   * @brief Configure whether input arrays are partial (for distributed processing)
   * @param b True if input arrays represent only a portion of the domain
   */
  void set_input_array_partial(bool b) {is_input_array_partial = b;}

  /**
   * @brief Enable/disable default domain partitioning
   * @param b True to use default domain partitioning
   */
  void set_use_default_domain_partition(bool b) {use_default_domain_partition = true;}

  /**
   * @brief Convert a tracker name string to tracker type ID
   * @param s Tracker name (e.g., "cp", "critical_point", "iso", "contour")
   * @return Tracker type ID, or 0 if name is not recognized
   */
  static int str2tracker(const std::string&);

public:
  /**
   * @brief Initialize the tracker
   *
   * Called once before processing any data. Sets up data structures,
   * domain decomposition, and any required preprocessing.
   */
  virtual void initialize() = 0;

  /**
   * @brief Finalize tracking and produce results
   *
   * Called after all timesteps have been processed. Completes trajectory
   * assembly and prepares final output data.
   */
  virtual void finalize() = 0;

  /**
   * @brief Advance to the next timestep
   * @return True if more timesteps remain, false if processing is complete
   */
  virtual bool advance_timestep() = 0;

  /**
   * @brief Process the current timestep
   *
   * Analyzes the current snapshot(s) to detect features and update
   * tracking state. Called once per timestep.
   */
  virtual void update_timestep() = 0;

public:
  /**
   * @brief Push a field data snapshot for the next timestep
   * @param g Shared pointer to ndarray_group containing field data
   */
  virtual void push_field_data_snapshot(std::shared_ptr<ndarray_group<>> g) {snapshots.push_back(g);}

  /**
   * @brief Push a single field data array as a snapshot
   * @param key Field name/key
   * @param arr Field data array
   */
  virtual void push_field_data_snapshot(const std::string key, const ndarray<double>& arr) {
    std::shared_ptr<ndarray_group<>> g(new ndarray_group<>);
    g->set(key, arr);
    push_field_data_snapshot(g);
  }

  /**
   * @brief Remove the oldest snapshot from the queue
   * @return True if a snapshot was removed, false if queue was empty
   */
  virtual bool pop_field_data_snapshot();

public:
  void set_fixed_quantization_factor(bool, double); // use 

protected:
  std::deque<std::shared_ptr<ndarray_group<>>> snapshots;
  int ntimesteps = 0; // unlimited

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
  else if (s == "mpas-o-pt" || s == "mpas-ocean-pt" || s == "mpas-ocean-particles")
    return TRACKER_MPAS_O_PARTICLES;
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
