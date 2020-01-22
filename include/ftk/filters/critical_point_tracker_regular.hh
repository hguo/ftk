#ifndef _FTK_CRITICAL_POINT_TRACKER_REGULAR_HH
#define _FTK_CRITICAL_POINT_TRACKER_REGULAR_HH

#include <ftk/ndarray.hh>
#include <ftk/hypermesh/lattice_partitioner.hh>
#include <ftk/filters/critical_point_tracker.hh>
#include <ftk/external/diy-ext/gather.hh>

namespace ftk {

enum {
  SOURCE_NONE, 
  SOURCE_GIVEN, // explicit
  SOURCE_DERIVED // implicit
};

struct critical_point_tracker_regular : public critical_point_tracker {
  critical_point_tracker_regular() {}
  critical_point_tracker_regular(int argc, char **argv) : critical_point_tracker(argc, argv) {}

  void set_domain(const lattice& l) {domain = l;} // spatial domain
  void set_array_domain(const lattice& l) {array_domain = l;}
  void set_start_timestep(int);
  void set_end_timestep(int);

  void set_use_default_domain_partition(bool b) {use_default_domain_partition = true;}
  void set_input_array_partial(bool b) {is_input_array_partial = b;}
  void set_local_domain(const lattice&); // rank-specific "core" region of the block
  void set_local_array_domain(const lattice&); // rank-specific "ext" region of the block

  void set_scalar_field_source(int s) {scalar_field_source = s;}
  void set_vector_field_source(int s) {vector_field_source = s;}
  void set_jacobian_field_source(int s) {jacobian_field_source = s;}
  void set_jacobian_symmetric(bool s) {is_jacobian_field_symmetric = s;}

  void set_type_filter(unsigned int);

  virtual void initialize() = 0;
  virtual void finalize() = 0;

  virtual void set_current_timestep(int t) {current_timestep = t;}
  virtual bool advance_timestep();
  virtual void update_timestep() = 0;

  void set_coordinates(const ndarray<double>& coords_) {coords = coords_; use_explicit_coords = true;}
  virtual void push_snapshot_scalar_field(const ndarray<double>& scalar0) {scalar.push_back(scalar0);}
  virtual void push_snapshot_vector_field(const ndarray<double>& V0) {V.push_back(V0);}
  virtual void push_snapshot_jacobian_field(const ndarray<double>& gradV0) {gradV.push_back(gradV0);}

  // pushing 3D/4D spacetime volume directly
  void push_spacetime_scalar_field(const ndarray<double>& scalar0);
  void push_spacetime_vector_field(const ndarray<double>& V0);
  void push_spacetime_jacobian_field(const ndarray<double>& gradV0);

  void set_max_num_staged_timesteps(int);

protected:
  template <int N, typename T=double>
  bool filter_critical_point_type(const critical_point_t<N, T>& cp);

protected: // config
  lattice domain, array_domain, 
          local_domain, local_array_domain;
  // lattice_partitioner partitioner;

  bool use_default_domain_partition = true;
  bool is_input_array_partial = false;

  int start_timestep = 0, 
      end_timestep = std::numeric_limits<int>::max();
  int max_num_staged_timesteps = 2;

  int scalar_field_source = SOURCE_NONE, 
      vector_field_source = SOURCE_NONE,
      jacobian_field_source = SOURCE_NONE;
  bool use_explicit_coords = false;
  bool is_jacobian_field_symmetric = false;
  bool use_type_filter = false;
  unsigned int type_filter = 0;

protected:
  ndarray<double> coords;
  std::deque<ndarray<double>> scalar, V, gradV;
  int current_timestep = 0;
};

/////
bool critical_point_tracker_regular::advance_timestep()
{
  update_timestep();

  // fprintf(stderr, "%lu, %lu\n", scalar.size(), V.size());
  scalar.pop_front();
  V.pop_front();
  gradV.pop_front();
#if 0
  if (scalar.size() > max_num_staged_timesteps) scalar.pop_back();
  if (V.size() > max_num_staged_timesteps) V.pop_back();
  if (gradV.size() > max_num_staged_timesteps) gradV.pop_back();
#endif

  current_timestep ++;
  return V.size() > 1; // > 0;
}

inline void critical_point_tracker_regular::set_type_filter(unsigned int f)
{
  use_type_filter = true;
  type_filter = f;
}
  
template <int N, typename T>
bool critical_point_tracker_regular::filter_critical_point_type(
    const critical_point_t<N, T>& cp)
{
  // fprintf(stderr, "typefilter=%lu, type=%lu\n", 
  //     type_filter, cp.type);
  if (use_type_filter) {
    if (type_filter & cp.type) return true;
    else return false;
  }
  else return true;
}

inline void critical_point_tracker_regular::push_spacetime_scalar_field(const ndarray<double>& scalar_all)
{
  for (size_t t = 0; t < scalar_all.shape(scalar_all.nd()-1); t ++) {
    auto array = scalar_all.slice_time(t);
    push_snapshot_scalar_field(array);
  }
}

}

#endif
