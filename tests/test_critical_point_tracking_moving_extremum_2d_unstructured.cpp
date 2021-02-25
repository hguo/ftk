#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include "constants.hh"
#include <ftk/filters/critical_point_tracker_2d_unstructured.hh>

using nlohmann::json;

const std::string mesh_filename = "1x1.vtu";
const int nt = 32;
const double kernel_size = 0.2;

#if FTK_HAVE_VTK
TEST_CASE("critical_point_tracking_moving_extremum_2d_unstructured") {
  diy::mpi::communicator world;
 
  ftk::simplicial_unstructured_2d_mesh<> m;
  m.from_vtk_unstructured_grid_file(mesh_filename);
  m.build_smoothing_kernel(kernel_size);

  ftk::critical_point_tracker_2d_unstructured tracker(world, m);
  tracker.initialize();

  for (int i = 0; i < nt; i ++) {
    auto data = ftk::synthetic_moving_extremum_unstructured<double, 2>(
        m.get_coords(), 
        {0.5, 0.5}, // center
        {0.02, 0.02}, // direction
        static_cast<double>(i) // time
    );

    ftk::ndarray<double> scalar, grad, J;
    m.smooth_scalar_gradient_jacobian(data, /*kernel_size,*/ scalar, grad, J);

    // write back for debugging
    char filename[1024];
    sprintf(filename, "moving_extremum-%02d.vtu", i);
    m.scalar_to_vtk_unstructured_grid_data_file(filename, "scalar", data);
    sprintf(filename, "moving_extremum-smooth-%02d.vtu", i);
    m.scalar_to_vtk_unstructured_grid_data_file(filename, "scalar", scalar);
    sprintf(filename, "moving_extremum_grad-%02d.vtu", i);
    m.vector_to_vtk_unstructured_grid_data_file(filename, "grad", grad);

    scalar.reshape({1, scalar.dim(0)});
    tracker.push_field_data_snapshot(scalar, grad, J);

    if (i != 0) tracker.advance_timestep();
    else tracker.update_timestep();
  }

  tracker.finalize();
  tracker.write_traced_critical_points_vtk("moving_extremum_2d_unstructured.vtp");

  auto trajs = tracker.get_traced_critical_points();
  
  if (world.rank() == 0)
    REQUIRE(trajs.size() == 1);
}
#endif

#include "main.hh"
