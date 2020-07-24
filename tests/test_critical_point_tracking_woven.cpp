#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include "constants.hh"
#include <ftk/filters/critical_point_tracker_wrapper.hh>

using nlohmann::json;

const int woven_n_trajs = 48;

TEST_CASE("critical_point_tracking_woven_synthetic") {
  auto result = track_cp2d(js_woven_synthetic);
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == woven_n_trajs);
    // REQUIRE(std::get<1>(result) == 4194); // TODO: check out if this number varies over different platform
}

TEST_CASE("critical_point_tracking_woven_float64") {
  auto result = track_cp2d(js_woven_float64);
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == woven_n_trajs);
}

#if FTK_HAVE_NETCDF
TEST_CASE("critical_point_tracking_woven_nc") {
  auto result = track_cp2d(js_woven_nc_unlimited_time);
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == woven_n_trajs);
}

TEST_CASE("critical_point_tracking_woven_nc_no_time") {
  auto result = track_cp2d(js_woven_nc_no_time);
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == woven_n_trajs);
}
#endif

#if FTK_HAVE_VTK
TEST_CASE("critical_point_tracking_woven_vti") {
  auto result = track_cp2d(js_woven_vti);
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == woven_n_trajs);
}
#endif

int main(int argc, char **argv)
{
  diy::mpi::environment env;
  
  // auto result = track_cp2d(js_woven_float64);
  // auto result = track_cp2d(js_woven_nc_no_time);
  // return 0;

  Catch::Session session;
  return session.run(argc, argv);
}
