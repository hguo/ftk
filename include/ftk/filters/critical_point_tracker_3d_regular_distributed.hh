#ifndef _FTK_CRITCAL_POINT_TRACKER_3D_REGULAR_DISTRIBUTED_HH
#define _FTK_CRITCAL_POINT_TRACKER_3D_REGULAR_DISTRIBUTED_HH

#include <ftk/filters/critical_point_tracker_3d_regular_streaming.hh>
#include <ftk/hypermesh/lattice_partitioner.hh>
#include <ftk/external/diy/master.hpp>
#include <ftk/external/diy-ext/gather.hh>
#include <mpi.h>

namespace ftk {

struct critical_point_tracker_3d_regular_distributed
  : public critical_point_tracker_3d_regular
{
  critical_point_tracker_3d_regular_distributed()
    : partitioner(m.get_lattice()) {}
  virtual ~critical_point_tracker_3d_regular_distributed() {};

  void update();

  void set_input_scalar_field_distributed(const ndarray<double>& scalar); // get the data block
  void set_input_vector_field_distributed(const ndarray<double>& V);
  void set_input_jacobian_field_distributed(const ndarray<double>& gradV);

#if 0
  void simplex_vectors(const std::vector<std::vector<int>>& vertices, double v[3][2]) const;
  void simplex_scalars(const std::vector<std::vector<int>>& vertices, double values[3]) const;
  void simplex_jacobians(const std::vector<std::vector<int>>& vertices, 
      double Js[3][2][2]) const;
#endif

protected:
  ndarray<double> scalar_part, V_part, gradV_part; // TODO: use part data instead of full data
  lattice_partitioner partitioner;

  typedef regular_simplex_mesh_element element_t;
  // std::map<std::string, critical_point_3dt_t> discrete_critical_points_map;
};

//////////////////
void critical_point_tracker_3d_regular_distributed::update()
{
  // initializing vector fields
  if (has_scalar_field) {
    if (!has_vector_field) {
      V = gradient3Dt(scalar);
      has_vector_field = true;
    }
    if (!has_jacobian_field) {
      gradV = jacobian3Dt(V);
      has_jacobian_field = true;
    }
    symmetric_jacobian = true;
  }

  // initializing bounds
  if (m.lb() == m.ub()) {
    if (has_scalar_field) // default lb/ub for scalar field
      m.set_lb_ub({2, 2, 2, 0}, {
          static_cast<int>(V.dim(1)-3), 
          static_cast<int>(V.dim(2)-3), 
          static_cast<int>(V.dim(3)-3),
          static_cast<int>(V.dim(4)-1)});
    else // defaulat lb/ub for vector field
      m.set_lb_ub({0, 0, 0, 0}, {
          static_cast<int>(V.dim(1)-1), 
          static_cast<int>(V.dim(2)-1), 
          static_cast<int>(V.dim(3)-1), 
          static_cast<int>(V.dim(4)-1)});
  }
 
  // initialize partitions
  partitioner.partition(comm.size(), 
      {}, // no given number of cuts
      {1, 1, 1, 0} // ghost size is 1 for spatial dimensions
  );
  // std::cerr << partitioner << std::endl;
  
  // scan 2-simplices
  if (comm.rank() == 0) fprintf(stderr, "tracking 3D critical points...\n");
  m.element_for(3,
      partitioner.get_ext(comm.rank()), 
      ftk::ELEMENT_SCOPE_ALL, 
      [=](element_t e) {
        critical_point_3dt_t cp;
        if (check_simplex(e, cp)) {
          std::lock_guard<std::mutex> guard(mutex);
          discrete_critical_points[e] = cp;
          // fprintf(stderr, "%f, %f, %f, %f, type=%d\n", cp[0], cp[1], cp[2], cp[3], cp.type);
        }
      });

  // gathering all discrete critical points
  diy::mpi::gather(comm, discrete_critical_points, discrete_critical_points, 0);

  if (comm.rank() == 0) {
    // fprintf(stderr, "trace intersections... %lu\n", discrete_critical_points.size());
    // trace_intersections();

    // convert connected components to traced critical points
    fprintf(stderr, "tracing critical points...\n");
    trace_connected_components();
  }
}

#if 0
void critical_point_tracker_3d_regular_distributed::simplex_vectors(const std::vector<std::vector<int>>& vertices, double v[3][2]) const
{
  // TODO
}

void critical_point_tracker_3d_regular_distributed::simplex_scalars(const std::vector<std::vector<int>>& vertices, double values[3]) const
{
  // TODO
}

void critical_point_tracker_3d_regular_distributed::simplex_jacobians(
    const std::vector<std::vector<int>>& vertices, 
    double Js[3][2][2]) const
{
  // TODO
}
#endif

}

#endif
