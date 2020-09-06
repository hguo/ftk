#ifndef _FTK_PARALLEL_VECTOR_TRACKER_HH
#define _FTK_PARALLEL_VECTOR_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/filters/filter.hh>
#include <ftk/filters/parallel_vector_curve_set.hh>

namespace ftk {

struct parallel_vector_tracker_3d_regular : public filter 
{
  parallel_vector_tracker_3d_regular();
  
  void set_domain(const lattice& l) {domain = l;} // spatial domain
  void set_array_domain(const lattice& l) {array_domain = l;}

  void push_field_data_snapshot(
      const ndarray<double> &v, 
      const ndarray<double> &w);

protected:
  simplicial_regular_mesh m;
  typedef simplicial_regular_mesh_element element_t;
  typedef parallel_vector_point_t pv_t;

  std::map<element_t, parallel_vector_point_t> discrete_pvs; // discrete parallel vector points
  
protected:
  lattice domain, array_domain;
};

}

#endif
