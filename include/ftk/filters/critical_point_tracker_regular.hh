#ifndef _FTK_CRITICAL_POINT_TRACKER_REGULAR_HH
#define _FTK_CRITICAL_POINT_TRACKER_REGULAR_HH

#include <ndarray/ndarray.hh>
#include <ndarray/lattice_partitioner.hh>
#include <ftk/filters/critical_point_tracker.hh>
#include <ftk/filters/regular_tracker.hh>
#include <ftk/utils/gather.hh>

namespace ftk {

/**
 * @brief Critical point tracker for regular (structured) grids
 *
 * This is an abstract base class that extends critical_point_tracker for
 * regular/structured grids. It combines the general critical point tracking
 * functionality with optimizations specific to regular grid topologies.
 *
 * Regular grids have implicit connectivity, which enables more efficient
 * algorithms for:
 * - Neighbor lookups (no explicit adjacency lists needed)
 * - Domain decomposition (regular block partitioning)
 * - GPU acceleration (regular memory access patterns)
 *
 * This class is not meant for direct use by end users. Use specific
 * implementations like critical_point_tracker_2d_regular or
 * critical_point_tracker_3d_regular instead.
 */
struct critical_point_tracker_regular : public critical_point_tracker, public regular_tracker {
  critical_point_tracker_regular(diy::mpi::communicator comm, int nd) : critical_point_tracker(comm), regular_tracker(comm, nd), tracker(comm) {}
  virtual ~critical_point_tracker_regular() {}

protected:
  typedef simplicial_regular_mesh_element element_t;

  std::map<element_t, feature_point_t> discrete_critical_points;
  std::vector<std::set<element_t>> connected_components;

public: // cp io
  /**
   * @brief Get discrete critical points (not yet assembled into trajectories)
   * @return Map from mesh elements to critical points
   */
  const std::map<element_t, feature_point_t>& get_discrete_critical_points() const {return discrete_critical_points;}

  /**
   * @brief Get critical points as a vector
   * @return Vector of all detected critical points
   */
  std::vector<feature_point_t> get_critical_points() const override;
  // void put_critical_points(const std::vector<feature_point_t>&);
};

/////
////
inline std::vector<feature_point_t> critical_point_tracker_regular::get_critical_points() const
{
  std::vector<feature_point_t> results;
  for (const auto &kv : discrete_critical_points) 
    results.push_back(kv.second);
  return results;
}

}

#endif
