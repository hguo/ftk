#ifndef _FTK_CONTOUR_TRACKER_REGULAR_HH
#define _FTK_CONTOUR_TRACKER_REGULAR_HH

#include <ftk/ndarray.hh>
#include <ftk/mesh/lattice_partitioner.hh>
#include <ftk/filters/contour_tracker.hh>
#include <ftk/filters/regular_tracker.hh>
#include <ftk/utils/gather.hh>

namespace ftk {

// this is an abstract class, not for users
struct contour_tracker_regular : public contour_tracker, public regular_tracker {
  contour_tracker_regular(diy::mpi::communicator comm, int nd/*2 or 3*/) : contour_tracker(comm), regular_tracker(comm, nd), tracker(comm) {}
  virtual ~contour_tracker_regular() {}

  void reset();

protected:
  typedef simplicial_regular_mesh_element element_t;
  
  std::map<element_t, feature_point_t> intersections;
  std::set<element_t> related_cells;

protected:
  void build_isovolumes();

public: // cp io
  const std::map<element_t, feature_point_t>& get_discrete_intersections() const {return intersections;}
  std::vector<feature_point_t> get_intersections() const;
};

inline void contour_tracker_regular::reset()
{
  current_timestep = 0;
  intersections.clear();
  related_cells.clear();
  contour_tracker::reset();
}

/////
////
// inline std::vector<contour_t> contour_tracker_regular::get_contours() const
// {
//   std::vector<contour_t> results;
//   for (const auto &kv : discrete_contours) 
//     results.push_back(kv.second);
//   return results;
// }

// inline void contour_tracker_regular::put_contours(const std::vector<contour_t>& data) 
// {
//   for (const auto& cp : data) {
//     element_t e(m, cpdims(), cp.tag);
//     discrete_contours[e] = cp;
//   }
// }

inline std::vector<feature_point_t> contour_tracker_regular::get_intersections() const
{
  const auto pts = get_discrete_intersections();
  std::vector<feature_point_t> results;

  for (const auto &kv : pts)
    results.push_back(kv.second);

  return results;
}

}

#endif
