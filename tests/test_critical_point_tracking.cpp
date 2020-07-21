#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include "test_constants.hh"
#include <ftk/filters/critical_point_tracker_wrapper.hh>

using nlohmann::json;

std::tuple<size_t, size_t> track2D(const json& jstream)
{
  ftk::ndarray_stream<> stream;
  stream.configure(jstream);

  ftk::critical_point_tracker_wrapper consumer;
  consumer.consume(stream);

  diy::mpi::communicator world;
  if (world.rank() == 0) {
    auto tracker = std::dynamic_pointer_cast<ftk::critical_point_tracker_2d_regular>( consumer.get_tracker() );
    auto trajs = tracker->get_traced_critical_points();
    auto points = tracker->get_discrete_critical_points();
    return {trajs.size(), points.size()};
  } else 
    return {0, 0};
}

TEST_CASE("critical_point_tracking_woven_synthetic") {
  auto result = track2D(js_woven_synthetic);
  REQUIRE(std::get<0>(result) == 48);
  // REQUIRE(std::get<1>(result) == 4194); // TODO: check out if this number varies over different platform
}

TEST_CASE("critical_point_tracking_woven_float64") {
  auto result = track2D(js_woven_float64);
  REQUIRE(std::get<0>(result) == 48);
}

#if FTK_HAVE_NETCDF
// TEST_CASE("critical_point_tracking_woven_nc") {
//   auto result = track2D(js_woven_nc);
//   REQUIRE(std::get<0>(result) == 48);
// }

TEST_CASE("critical_point_tracking_woven_nc_no_time") {
  auto result = track2D(js_woven_nc_no_time);
  REQUIRE(std::get<0>(result) == 48);
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
