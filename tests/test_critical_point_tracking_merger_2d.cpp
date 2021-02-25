#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include "constants.hh"

#include "main.hh"

using nlohmann::json;

const int merger_n_trajs = 9;
  
TEST_CASE("critical_point_tracking_merger_2d_write_discrete_critical_points_binary") {
  auto result = track_cp2d(js_merger_2d_synthetic, {
    {"output", "merger_2d_discrete.bin"},
    {"output_type", "discrete"}
  });

  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == merger_n_trajs);
}

TEST_CASE("critical_point_tracking_merger_2d_read_discrete_critical_points_binary") {
  auto result = track_cp2d(js_merger_2d_synthetic, {
    {"archived_discrete_critical_points_filename", "merger_2d_discrete.bin"}
  });
  
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == merger_n_trajs);
}

TEST_CASE("critical_point_tracking_merger_2d_write_discrete_critical_points_json") {
  auto result = track_cp2d(js_merger_2d_synthetic, {
    {"output", "merger_2d_discrete.json"},
    {"output_type", "discrete"}
  });

  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == merger_n_trajs);
}

TEST_CASE("critical_point_tracking_merger_2d_read_discrete_critical_points_json") {
  auto result = track_cp2d(js_merger_2d_synthetic, {
    {"archived_discrete_critical_points_filename", "merger_2d_discrete.json"}
  });
  
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == merger_n_trajs);
}

/////
TEST_CASE("critical_point_tracking_merger_2d_write_traced_critical_points_binary") {
  auto result = track_cp2d(js_merger_2d_synthetic, {
    {"output", "merger_2d_traced.bin"},
    {"output_type", "traced"}
  });

  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == merger_n_trajs);
}

TEST_CASE("critical_point_tracking_merger_2d_read_traced_critical_points_binary") {
  auto result = track_cp2d(js_merger_2d_synthetic, {
    {"archived_traced_critical_points_filename", "merger_2d_traced.bin"}
  });
  
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == merger_n_trajs);
}

TEST_CASE("critical_point_tracking_merger_2d_write_traced_critical_points_json") {
  auto result = track_cp2d(js_merger_2d_synthetic, {
    {"output", "merger_2d_traced.json"},
    {"output_type", "traced"}
  });

  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == merger_n_trajs);
}

TEST_CASE("critical_point_tracking_merger_2d_read_traced_critical_points_json") {
  auto result = track_cp2d(js_merger_2d_synthetic, {
    {"archived_traced_critical_points_filename", "merger_2d_traced.json"}
  });
  
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == merger_n_trajs);
}

