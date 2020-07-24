#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include "constants.hh"
#include <ftk/filters/critical_point_tracker_wrapper.hh>

using nlohmann::json;

TEST_CASE("critical_point_tracking_moving_extremum_2d") {
  auto result = track_cp2d(js_moving_extremum_2d_vti);
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

    if (world.rank() == 0) {
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
}

int main(int argc, char **argv)
{
  diy::mpi::environment env;
  
  Catch::Session session;
  return session.run(argc, argv);
}
