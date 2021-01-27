#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <ftk/mesh/simplicial_unstructured_extruded_3d_mesh.hh>
#include <ftk/mesh/simplicial_xgc_3d_mesh.hh>
#include <ftk/ndarray.hh>

using namespace ftk;

diy::mpi::environment env;

int nphi = 16, iphi = 1;
const int vphi = 8;

std::string xgc_data_path; 
std::string xgc_mesh_filename = "xgc.mesh.h5", 
            xgc_bfield_filename = "xgc.bfield.h5",
            xgc_3d_filename = "xgc.3d.00001.h5",
            xgc_units_filename = "units.m";

const double xgc_smoothing_kernel_size = 0.03;

std::shared_ptr<simplicial_xgc_2d_mesh<>> mx2; // 2d mesh
std::shared_ptr<simplicial_xgc_3d_mesh<>> mx3, // 3d mesh, 
                                          mx30; // 3d mesh for write-backs

#if FTK_HAVE_HDF5
TEST_CASE("xgc_synthetic_filament_tracking") {
  mx2 = simplicial_xgc_2d_mesh<>::from_xgc_mesh_h5(xgc_data_path + "/" + xgc_mesh_filename);
  mx2->initialize_point_locator();
  mx2->read_bfield_h5(xgc_data_path + "/" + xgc_bfield_filename);
  mx2->read_units_m(xgc_data_path + "/" + xgc_units_filename);
  mx2->build_smoothing_kernel(xgc_smoothing_kernel_size);
  
  mx3.reset( new simplicial_xgc_3d_mesh<>(mx2) );
  mx3->probe_nphi_iphi_h5( xgc_data_path + "/" + xgc_3d_filename );
  mx3->set_vphi( vphi );
  mx3->initialize_interpolants();
}
#endif

int main(int argc, char **argv)
{
  const char* path = std::getenv("FTK_XGC_TEST_DATA_PATH");
  if (path)
    xgc_data_path = path;

  Catch::Session session;
  session.run(argc, argv);
  return 0;
}
