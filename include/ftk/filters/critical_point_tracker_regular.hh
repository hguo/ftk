#ifndef _FTK_CRITICAL_POINT_TRACKER_REGULAR_HH
#define _FTK_CRITICAL_POINT_TRACKER_REGULAR_HH

#include <ftk/ndarray.hh>
#include <ftk/mesh/lattice_partitioner.hh>
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
  virtual ~critical_point_tracker_regular() {}

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

  virtual void push_scalar_field_snapshot(const ndarray<double>&) = 0;
  virtual void push_vector_field_snapshot(const ndarray<double>&) = 0;

  virtual void initialize() = 0;
  virtual void finalize() = 0;

  virtual bool advance_timestep();
  virtual void update_timestep() = 0;

  void set_coordinates(const ndarray<double>& coords_) {coords = coords_; use_explicit_coords = true;}

protected:
  bool filter_critical_point_type(const critical_point_t& cp);

protected: // config
  lattice domain, array_domain, 
          local_domain, local_array_domain;
  // lattice_partitioner partitioner;

  bool use_default_domain_partition = true;
  bool is_input_array_partial = false;

  int start_timestep = 0, 
      end_timestep = std::numeric_limits<int>::max();

  int scalar_field_source = SOURCE_NONE, 
      vector_field_source = SOURCE_NONE,
      jacobian_field_source = SOURCE_NONE;
  bool use_explicit_coords = false;
  bool is_jacobian_field_symmetric = false;

protected:
  ndarray<double> coords;
  // std::deque<ndarray<double>> scalar, V, gradV;
};

/////
inline bool critical_point_tracker_regular::advance_timestep()
{
  update_timestep();
  pop_field_data_snapshot();

  current_timestep ++;
  return field_data_snapshots.size() > 0;
}
  
inline bool critical_point_tracker_regular::filter_critical_point_type(
    const critical_point_t& cp)
{
  // fprintf(stderr, "typefilter=%lu, type=%lu\n", 
  //     type_filter, cp.type);
  if (use_type_filter) {
    if (type_filter & cp.type) return true;
    else return false;
  }
  else return true;
}

}

#endif
