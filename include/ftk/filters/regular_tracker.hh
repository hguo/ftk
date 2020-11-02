#ifndef _FTK_REGULAR_TRACKER_HH
#define _FTK_REGULAR_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/mesh/lattice_partitioner.hh>
#include <ftk/filters/filter.hh>

namespace ftk {

struct regular_tracker : public virtual filter {
  regular_tracker(int nd/*2 or 3*/) : m(nd+1) {}
  virtual ~regular_tracker() {}
 
public:
  void set_domain(const lattice& l) {domain = l;} // spatial domain
  void set_array_domain(const lattice& l) {array_domain = l;}

  void set_local_domain(const lattice&); // rank-specific "core" region of the block
  void set_local_array_domain(const lattice&); // rank-specific "ext" region of the block

  void set_coordinates(const ndarray<double>& coords_) {coords = coords_; use_explicit_coords = true;}

protected:
  simplicial_regular_mesh m;
  
  lattice domain, array_domain, 
          local_domain, local_array_domain;

  bool use_explicit_coords = false;
  ndarray<double> coords;

protected: // internal use
  template <typename I=int> void simplex_indices(const std::vector<std::vector<int>>& vertices, I indices[]) const;
};

/////////////////////////////
template <typename I>
inline void regular_tracker::simplex_indices(
    const std::vector<std::vector<int>>& vertices, I indices[]) const
{
  for (int i = 0; i < vertices.size(); i ++)
    indices[i] = m.get_lattice().to_integer(vertices[i]);
}

}


#endif
