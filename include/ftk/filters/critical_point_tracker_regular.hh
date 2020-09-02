#ifndef _FTK_CRITICAL_POINT_TRACKER_REGULAR_HH
#define _FTK_CRITICAL_POINT_TRACKER_REGULAR_HH

#include <ftk/ndarray.hh>
#include <ftk/mesh/lattice_partitioner.hh>
#include <ftk/filters/critical_point_tracker.hh>
#include <ftk/external/diy-ext/gather.hh>

namespace ftk {

// this is an abstract class, not for users
struct critical_point_tracker_regular : public critical_point_tracker {
  critical_point_tracker_regular(int nd/*2 or 3*/) : m(nd+1) {}
  virtual ~critical_point_tracker_regular() {}

  void set_domain(const lattice& l) {domain = l;} // spatial domain
  void set_array_domain(const lattice& l) {array_domain = l;}

  void set_local_domain(const lattice&); // rank-specific "core" region of the block
  void set_local_array_domain(const lattice&); // rank-specific "ext" region of the block

  void set_coordinates(const ndarray<double>& coords_) {coords = coords_; use_explicit_coords = true;}

protected:
  simplicial_regular_mesh m;
  typedef simplicial_regular_mesh_element element_t;
  
  std::map<element_t, critical_point_t> discrete_critical_points;
  std::vector<std::set<element_t>> connected_components;

public: // cp io
  const std::map<element_t, critical_point_t>& get_discrete_critical_points() const {return discrete_critical_points;}

  std::vector<critical_point_t> get_critical_points() const;
  void put_critical_points(const std::vector<critical_point_t>&);

protected: // internal use
  template <typename I=int> void simplex_indices(const std::vector<std::vector<int>>& vertices, I indices[]) const;

protected: // config
  lattice domain, array_domain, 
          local_domain, local_array_domain;
  // lattice_partitioner partitioner;

  bool use_explicit_coords = false;

protected:
  ndarray<double> coords;
  // std::deque<ndarray<double>> scalar, V, gradV;
};

/////
////
inline std::vector<critical_point_t> critical_point_tracker_regular::get_critical_points() const
{
  std::vector<critical_point_t> results;
  for (const auto &kv : discrete_critical_points) 
    results.push_back(kv.second);
  return results;
}

inline void critical_point_tracker_regular::put_critical_points(const std::vector<critical_point_t>& data) 
{
  for (const auto& cp : data) {
    element_t e(m, cpdims(), cp.tag);
    discrete_critical_points[e] = cp;
  }
}

template <typename I>
inline void critical_point_tracker_regular::simplex_indices(
    const std::vector<std::vector<int>>& vertices, I indices[]) const
{
  for (int i = 0; i < vertices.size(); i ++)
    indices[i] = m.get_lattice().to_integer(vertices[i]);
}

}

#endif
