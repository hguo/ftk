#ifndef _FTK_CRITICAL_POINT_TRACKER_2D_REGULAR_DISTRIBUTED_STREAMING_HH
#define _FTK_CRITICAL_POINT_TRACKER_2D_REGULAR_DISTRIBUTED_STREAMING_HH

#include <ftk/filters/critical_point_tracker_2d_regular_distributed.hh>

namespace ftk {

struct critical_point_tracker_2d_regular_distributed_streaming
  : public critical_point_tracker_2d_regular_streaming
{
  void update_timestep();
  void update();
};

/////
void critical_point_tracker_2d_regular_distributed_streaming::update()
{
  std::vector<int> lb, ub;
  m.get_lb_ub(lb, ub);
  ub[2] = current_timestep - 1;
  m.set_lb_ub(lb, ub);
  
  // gathering all discrete critical points
  diy::mpi::gather(comm, discrete_critical_points, discrete_critical_points, 0);

  if (comm.rank() == 0) {
    // fprintf(stderr, "trace1\n");
    // trace_intersections();
    fprintf(stderr, "trace2\n");
    trace_connected_components();
    fprintf(stderr, "done.\n");
  }
}

void critical_point_tracker_2d_regular_distributed_streaming::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "processing timestep %d..\n", current_timestep);

  // initializing vector fields
  if (has_scalar_field) { 
    if (!has_vector_field) push_input_vector_field(gradient2D(scalar[0])); // 0 is the current timestep; 1 is the last timestep
    if (!has_jacobian_field) push_input_jacobian_field(jacobian2D(V[0]));
    symmetric_jacobian = true;
  }

  // initializing bounds
  if (m.lb() == m.ub()) {
    if (!scalar.empty())
      m.set_lb_ub({2, 2, 0}, {
          static_cast<int>(V[0].dim(1)-3), 
          static_cast<int>(V[0].dim(2)-3), 
          std::numeric_limits<int>::max()});
    else
      m.set_lb_ub({0, 0, 0}, {
          static_cast<int>(V[0].dim(1)-1), 
          static_cast<int>(V[0].dim(2)-1), 
          std::numeric_limits<int>::max()});
  }

  auto starts = m.get_lattice().starts(), 
       sizes = m.get_lattice().sizes();
  starts[2] = current_timestep;
  sizes[2] = 1;
  lattice current_lattice(starts, sizes);

  lattice_partitioner partitioner(current_lattice);
  partitioner.partition(comm.size(), {0, 0, 1}, {1, 1, 0});
  std::cerr << partitioner << std::endl;

  // scan 2-simplices
  // fprintf(stderr, "tracking 2D critical points...\n");
  auto func2 = [=](element_t e) {
      critical_point_2dt_t cp;
      if (check_simplex(e, cp)) {
        std::lock_guard<std::mutex> guard(mutex);
        fprintf(stderr, "%f, %f, %f\n", cp[0], cp[1], cp[2]);
        discrete_critical_points[e] = cp;
      }
    };

  if (V.size() >= 2) {
    // m.element_for_interval(2, current_timestep-1, current_timestep, func2, nthreads);
    auto local_lattice = partitioner.get_ext(comm.rank());
    local_lattice.starts_[2] = current_timestep - 1;
    m.element_for(2, local_lattice,
        ftk::ELEMENT_SCOPE_INTERVAL, 
        func2, nthreads);
  }

  // m.element_for_ordinal(2, current_timestep, func2, nthreads);
  m.element_for(2, 
      partitioner.get_ext(comm.rank()), 
      ftk::ELEMENT_SCOPE_ORDINAL,
      func2, nthreads);
}

}

#endif
