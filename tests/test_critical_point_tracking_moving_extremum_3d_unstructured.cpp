#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include "constants.hh"
#include <ftk/filters/critical_point_tracker_3d_unstructured.hh>

using nlohmann::json;
  
const std::string mesh_filename = "3d.vtu";
const int nt = 32;
const double kernel_size = 2.0;

#if FTK_HAVE_VTK
TEST_CASE("critical_point_tracking_moving_extremum_3d_unstructured") {
  diy::mpi::communicator world;
 
  ftk::simplicial_unstructured_3d_mesh<> m;
  m.from_vtk_unstructured_grid_file(mesh_filename);
  // m.build_smoothing_kernel(kernel_size);

  ftk::critical_point_tracker_3d_unstructured tracker(world, m);
  // tracker.set_number_of_threads(1);
  tracker.initialize();

  for (int i = 0; i < nt; i ++) {
    ftk::ndarray<double> grad = ftk::synthetic_moving_extremum_grad_unstructured<double, 3>(
        m.get_coords(), 
        {0.0, 0.0, 0.0}, // center
        // {0.1, 0.1, 0.1}, // direction
        {1.0, 1.0, 1.0}, // direction
        // {-0.1, 0.1, -0.2}, // direction
        static_cast<double>(i) // time
    );
   
    // write back
    char filename[1024];
    sprintf(filename, "moving_extremum-3d-grad-%03d.vtu", i);
    m.vector_to_vtk_unstructured_grid_data_file(filename, "grad", grad);
    
    tracker.push_field_data_snapshot(ftk::ndarray<double>(), grad, ftk::ndarray<double>());

#if 0
    auto data = ftk::synthetic_moving_extremum_unstructured<double, 3>(
        m.get_coords(), 
        {0.0, 0.0, 0.0}, // center
        {1.0, 1.0, 1.0}, // direction
        static_cast<double>(i) // time
    );

    ftk::ndarray<double> scalar, grad, J;
    m.smooth_scalar_gradient_jacobian(data, kernel_size, scalar, grad, J);

    // write back for debugging
    char filename[1024];
    sprintf(filename, "moving_extremum-3d-%03d.vtu", i);
    m.scalar_to_vtk_unstructured_grid_data_file(filename, "scalar", data);
    sprintf(filename, "moving_extremum-3d-smooth-%03d.vtu", i);
    m.scalar_to_vtk_unstructured_grid_data_file(filename, "scalar", scalar);
    sprintf(filename, "moving_extremum-3d-grad-%03d.vtu", i);
    m.vector_to_vtk_unstructured_grid_data_file(filename, "grad", grad);

    scalar.reshape({1, scalar.dim(0)});
    tracker.push_field_data_snapshot(scalar, grad, J);
#endif

    if (i != 0) tracker.advance_timestep();
    else tracker.update_timestep();
  }

  tracker.finalize();
  tracker.write_traced_critical_points_vtk("moving_extremum_3d_unstructured.vtp");

  auto trajs = tracker.get_traced_critical_points();
  
  if (world.rank() == 0)
    REQUIRE(trajs.size() == 1);
}
#endif

#include "main.hh"
