#ifndef _FTK_CRITICAL_POINT_TRACKER_REGULAR_HH
#define _FTK_CRITICAL_POINT_TRACKER_REGULAR_HH

#include <ftk/ndarray.hh>
#include <ftk/mesh/lattice_partitioner.hh>
#include <ftk/filters/contour_tracker.hh>
#include <ftk/external/diy-ext/gather.hh>

namespace ftk {

// this is an abstract class, not for users
struct contour_tracker_regular : public contour_tracker {
  contour_tracker_regular(int nd/*2 or 3*/) : m(nd+1) {}
  virtual ~contour_tracker_regular() {}

  void set_domain(const lattice& l) {domain = l;} // spatial domain
  void set_array_domain(const lattice& l) {array_domain = l;}

protected:
  simplicial_regular_mesh m;
  typedef simplicial_regular_mesh_element element_t;
  
  std::map<element_t, feature_point_t> intersections;
  std::set<element_t> related_cells;

public: // cp io
  const std::map<element_t, feature_point_t>& get_discrete_intersections() const {return intersections;}

  std::vector<feature_point_t> get_intersections() const;

protected: // internal use
  template <typename I=int> void simplex_indices(const std::vector<std::vector<int>>& vertices, I indices[]) const;

protected: // config
  lattice domain, array_domain;
};

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

template <typename I>
inline void contour_tracker_regular::simplex_indices(
    const std::vector<std::vector<int>>& vertices, I indices[]) const
{
  for (int i = 0; i < vertices.size(); i ++)
    indices[i] = m.get_lattice().to_integer(vertices[i]);
}

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
