#ifndef _FTK_CRITCAL_POINT_TRACKER_2D_REGULAR_DISTRIBUTED_HH
#define _FTK_CRITCAL_POINT_TRACKER_2D_REGULAR_DISTRIBUTED_HH

#include <ftk/filters/critical_point_tracker_2d_regular_streaming.hh>
#include <ftk/hypermesh/lattice_partitioner.hh>
#include <ftk/external/diy/master.hpp>

namespace ftk {

struct critical_point_tracker_2d_regular_distributed
  : public critical_point_tracker_2d_regular
{
  critical_point_tracker_2d_regular_distributed()
    : partitioner(m.get_lattice()) {}
  virtual ~critical_point_tracker_2d_regular_distributed() {};

  void update();

  void set_input_scalar_field_distributed(const ndarray<double>& scalar); // get the data block
  void set_input_vector_field_distributed(const ndarray<double>& V);
  void set_input_jacobian_field_distributed(const ndarray<double>& gradV);

  void simplex_vectors(const std::vector<std::vector<int>>& vertices, double v[3][2]) const;
  void simplex_scalars(const std::vector<std::vector<int>>& vertices, double values[3]) const;
  void simplex_jacobians(const std::vector<std::vector<int>>& vertices, 
      double Js[3][2][2]) const;

protected:
  ndarray<double> scalar_part, V_part, gradV_part;
  lattice_partitioner partitioner;
  
  typedef regular_simplex_mesh_element element_t;
};
 
void critical_point_tracker_2d_regular_distributed::update()
{

}

void critical_point_tracker_2d_regular_distributed::simplex_vectors(const std::vector<std::vector<int>>& vertices, double v[3][2]) const
{
  // TODO
}

void critical_point_tracker_2d_regular_distributed::simplex_scalars(const std::vector<std::vector<int>>& vertices, double values[3]) const
{
  // TODO
}

void critical_point_tracker_2d_regular_distributed::simplex_jacobians(
    const std::vector<std::vector<int>>& vertices, 
    double Js[3][2][2]) const
{
  // TODO
}

}

#endif
