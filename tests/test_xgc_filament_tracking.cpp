#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <ftk/mesh/simplicial_unstructured_extruded_3d_mesh.hh>
#include <ftk/mesh/simplicial_xgc_3d_mesh.hh>
#include <ftk/filters/xgc_blob_filament_tracker.hh>
#include "main.hh"

using namespace ftk;

int nphi = 16, iphi = 1;
const int vphi = 8;

const std::string xgc_mesh_filename = "xgc.mesh.h5", 
      xgc_bfield_filename = "xgc.bfield.h5",
      xgc_3d_filename = "xgc.3d.00001.h5",
      xgc_units_filename = "units.m",
      xgc_synthetic_filename = "xgc.synthetic.*.h5",
      xgc_varname = "/dnOvernXGC";

const double xgc_smoothing_kernel_size = 0.03;

std::shared_ptr<simplicial_xgc_2d_mesh<>> mx2; // 2d mesh
std::shared_ptr<simplicial_xgc_3d_mesh<>> mx3, // 3d mesh, 
                                          mx30; // 3d mesh for write-backs

#if FTK_HAVE_HDF5
TEST_CASE("xgc_synthetic_filament_tracking") {
  diy::mpi::communicator world;

  mx2 = simplicial_xgc_2d_mesh<>::from_xgc_mesh_file(xgc_data_path + "/" + xgc_mesh_filename);
  mx2->initialize_point_locator();
  mx2->read_bfield(xgc_data_path + "/" + xgc_bfield_filename);
  mx2->read_units_m(xgc_data_path + "/" + xgc_units_filename);
  mx2->build_smoothing_kernel(xgc_smoothing_kernel_size);
  
  mx3.reset( new simplicial_xgc_3d_mesh<>(mx2) );
  mx3->probe_nphi_iphi( xgc_data_path + "/" + xgc_3d_filename );
  mx3->set_vphi( vphi );
  mx3->initialize_interpolants();
  
  std::vector<std::string> filenames = ftk::glob(xgc_synthetic_filename);
  REQUIRE( filenames.size() == 5 ); // make sure synthetic data are generated correctly
  REQUIRE( mx3->get_nphi() == 16);
  REQUIRE( mx3->get_iphi() == 1);
  
  std::shared_ptr<xgc_blob_filament_tracker> tracker;
  tracker.reset(new xgc_blob_filament_tracker(world, mx3));
  
  tracker->set_end_timestep(filenames.size() - 1);
  // tracker->use_accelerator(accelerator);
  // tracker->set_number_of_threads(nthreads);
  tracker->initialize();

  for (int k = 0; k < filenames.size(); k ++) {
    auto data = ndarray<double>::from_h5(filenames[k], xgc_varname)
      .get_transpose();

    REQUIRE(data.dim(0) == 56980);
    REQUIRE(data.dim(1) == 16);

    tracker->push_field_data_snapshot( data );

    if (k != 0) tracker->advance_timestep();
    if (k == filenames.size() - 1) tracker->update_timestep();
  }
  tracker->finalize();
  
  const feature_surface_t &surface = tracker->get_surfaces();
  REQUIRE(surface.pts.size() == 20667);
  REQUIRE(surface.tris.size() == 44420);
}
#endif
