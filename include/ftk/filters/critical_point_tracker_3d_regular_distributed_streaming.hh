#ifndef _FTK_CRITICAL_POINT_TRACKER_3D_REGULAR_DISTRIBUTED_STREAMING_HH
#define _FTK_CRITICAL_POINT_TRACKER_3D_REGULAR_DISTRIBUTED_STREAMING_HH

#include <ftk/filters/critical_point_tracker_3d_regular_distributed.hh>

namespace ftk {

struct critical_point_tracker_3d_regular_distributed_streaming
  : public critical_point_tracker_3d_regular_streaming
{
  void update_timestep();
  void update();
};

/////
void critical_point_tracker_3d_regular_distributed_streaming::update()
{
  // gathering all discrete critical points to the root process
  diy::mpi::gather(comm, discrete_critical_points, discrete_critical_points, 0);

  if (comm.rank() == 0) {
    trace_connected_components();
    fprintf(stderr, "done.\n");
  }
}

void critical_point_tracker_3d_regular_distributed_streaming::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "processing timestep %d..\n", current_timestep);

  // initializing vector fields
  if (has_scalar_field) { 
    if (!has_vector_field) push_input_vector_field(gradient3D(scalar[0])); // 0 is the current timestep; 1 is the last timestep
    if (!has_jacobian_field) push_input_jacobian_field(jacobian3D(V[0]));
    symmetric_jacobian = true;
  }

  // initializing bounds
  if (m.lb() == m.ub()) {
    if (!scalar.empty())
      m.set_lb_ub({2, 2, 0}, {
          static_cast<int>(V[0].dim(1)-3), 
          static_cast<int>(V[0].dim(2)-3), 
          static_cast<int>(V[0].dim(3)-3), 
          std::numeric_limits<int>::max()});
    else
      m.set_lb_ub({0, 0, 0}, {
          static_cast<int>(V[0].dim(1)-1), 
          static_cast<int>(V[0].dim(2)-1), 
          static_cast<int>(V[0].dim(3)-1), 
          std::numeric_limits<int>::max()});
  }

  auto starts = m.get_lattice().starts(), 
       sizes = m.get_lattice().sizes();
  starts[3] = current_timestep;
  sizes[3] = 1;
  lattice current_lattice(starts, sizes);

  lattice_partitioner partitioner(current_lattice);
  // partitioner.partition(comm.size(), {0, 0, 1}, {1, 1, 0});
  partitioner.partition(comm.size(), {0, 0, 0, 1}, {0, 0, 0, 0});
  // std::cerr << partitioner << std::endl;

  // scan 3-simplices
  // fprintf(stderr, "tracking 3D critical points...\n");
  auto func3 = [=](element_t e) {
      critical_point_3dt_t cp;
      if (check_simplex(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        // fprintf(stderr, "%f, %f, %f\n", cp[0], cp[1], cp[2]);
        discrete_critical_points[e] = cp;
      }
    };

  if (V.size() >= 2) {
    // m.element_for_interval(2, current_timestep-1, current_timestep, func3, nthreads);
    auto local_lattice = partitioner.get_ext(comm.rank());
    local_lattice.starts_[3] = current_timestep - 1;
    m.element_for(3, local_lattice,
        ftk::ELEMENT_SCOPE_INTERVAL, 
        func3, nthreads);
  }

  // m.element_for_ordinal(2, current_timestep, func3, nthreads);
  m.element_for(3,
      partitioner.get_ext(comm.rank()), 
      ftk::ELEMENT_SCOPE_ORDINAL,
      func3, nthreads);
}

}

#endif
