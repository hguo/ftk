#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include "constants.hh"

using nlohmann::json;
  
#if FTK_HAVE_VTK
TEST_CASE("critical_point_tracking_double_gyre_vti") {
  auto result = track_cp2d(js_double_gyre_vti, {
    {"output", "double-gyre-%04d.vtp"},
    {"output_type", "sliced"}
  });
  
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == 2);
}
#endif

TEST_CASE("critical_point_tracking_double_gyre") {
  auto result = track_cp2d(js_double_gyre_synthetic, {
    {"output", "double-gyre.vtp"},
    {"output_type", "traced"}
  });
  
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == 2);
}

#include "main.hh"
