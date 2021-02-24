#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include "constants.hh"
#include "main.hh"

using nlohmann::json;

const int woven_n_trajs = 56; // 48;

#if FTK_TEST_CUDA
TEST_CASE("critical_point_tracking_cuda_woven_synthetic") {
  auto result = track_cp2d(js_woven_synthetic, {
    {"accelerator", "cuda"}
  });
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == 48); // the domain is slightly different from the CPU version.  will fix it in the future.
}
#endif

TEST_CASE("critical_point_tracking_woven_synthetic") {
  auto result = track_cp2d(js_woven_synthetic);
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == woven_n_trajs);
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
