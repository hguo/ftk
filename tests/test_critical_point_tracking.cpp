#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include "constants.hh"
#include <ftk/filters/critical_point_tracker_wrapper.hh>

using nlohmann::json;

const int woven_n_trajs = 48;

std::tuple<size_t, size_t> track2D(const json& jstream)
{
  ftk::ndarray_stream<> stream;
  stream.configure(jstream);

  ftk::critical_point_tracker_wrapper consumer;
  consumer.consume(stream);
    
  auto tracker = std::dynamic_pointer_cast<ftk::critical_point_tracker_2d_regular>( consumer.get_tracker() );
  auto trajs = tracker->get_traced_critical_points();
  auto points = tracker->get_discrete_critical_points();
  return {trajs.size(), points.size()};
}

TEST_CASE("critical_point_tracking_woven_synthetic") {
  auto result = track2D(js_woven_synthetic);
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == woven_n_trajs);
    // REQUIRE(std::get<1>(result) == 4194); // TODO: check out if this number varies over different platform
}

TEST_CASE("critical_point_tracking_woven_float64") {
  auto result = track2D(js_woven_float64);
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == woven_n_trajs);
}

#if FTK_HAVE_NETCDF
TEST_CASE("critical_point_tracking_woven_nc") {
  auto result = track2D(js_woven_nc_unlimited_time);
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == woven_n_trajs);
}

TEST_CASE("critical_point_tracking_woven_nc_no_time") {
  auto result = track2D(js_woven_nc_no_time);
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == woven_n_trajs);
}
#endif

#if FTK_HAVE_VTK
TEST_CASE("critical_point_tracking_woven_vti") {
  auto result = track2D(js_woven_vti);
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == woven_n_trajs);
}

TEST_CASE("critical_point_tracking_moving_extremum_2d") {
  auto result = track2D(js_moving_extremum_2d_vti);
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == 1);
}

TEST_CASE("critical_point_tracking_moving_extremum_2d_random_motion") {
  const int ncases = 10;
  const int width = 21, height = 21;
  const int x0 = 10, y0 = 10;
  diy::mpi::communicator world;
  
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> d{0, 1};

  for (int k = 0; k < ncases; k ++) {
    const double dir_x = d(gen), dir_y = d(gen);
    json js = js_moving_extremum_2d_synthetic;
    js["width"] = width;
    js["height"] = height;
    js["x0"] = {x0, y0};
    js["dir"] = {dir_x, dir_y};

    ftk::ndarray_stream<> stream;
    stream.configure(js);

    ftk::critical_point_tracker_wrapper consumer;
    consumer.consume(stream);

    auto tracker = std::dynamic_pointer_cast<ftk::critical_point_tracker_2d_regular>( consumer.get_tracker() );
    auto trajs = tracker->get_traced_critical_points();

    REQUIRE(trajs.size() == 1);
    
    for (auto i = 0; i < trajs[0].size(); i ++) {
      const auto &p = trajs[0][i];
      double x = x0 + dir_x * p[2], 
             y = y0 + dir_y * p[2];

      REQUIRE(p[0] == Approx(x));
      REQUIRE(p[1] == Approx(y));
    }
  }
}
#endif

int main(int argc, char **argv)
{
  diy::mpi::environment env;
  
  // auto result = track2D(js_woven_float64);
  // auto result = track2D(js_woven_nc_no_time);
  // return 0;

  Catch::Session session;
  return session.run(argc, argv);
}
