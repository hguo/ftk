#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include "constants.hh"
#include <ftk/filters/critical_point_tracker_2d_unstructured.hh>

using nlohmann::json;

const std::string mesh_filename = "1x1.vtu";
const int nt = 32;
const double kernel_size = 0.045;

#if FTK_HAVE_VTK
TEST_CASE("critical_point_tracking_woven_unstructured") {
  diy::mpi::communicator world;
 
  ftk::simplicial_unstructured_2d_mesh<> m;
  m.from_vtk_unstructured_grid_file(mesh_filename);
  m.build_smoothing_kernel(kernel_size);

  ftk::critical_point_tracker_2d_unstructured tracker(world, m);
  tracker.initialize();

  for (int i = 0; i < nt; i ++) {
    auto data = ftk::synthetic_woven_2D_unstructured<double>(m.get_coords(), i*0.1, 
        {0.5, 0.5}, // center
        10.0 // scaling factor
    );

    ftk::ndarray<double> scalar, grad, J;
    m.smooth_scalar_gradient_jacobian(data, scalar, grad, J);
    scalar.reshape({1, scalar.dim(0)});

    tracker.push_field_data_snapshot(scalar, grad, J);

    if (i != 0) tracker.advance_timestep();
    else tracker.update_timestep();
  }

  tracker.finalize();
  tracker.write_traced_critical_points_vtk("woven_unstructured.vtp");

  auto trajs = tracker.get_traced_critical_points();
  
  if (world.rank() == 0)
    REQUIRE(trajs.size() == 67); 
}
#endif

#include "main.hh"
