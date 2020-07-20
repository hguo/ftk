#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <ftk/filters/critical_point_tracker_wrapper.hh>

using nlohmann::json;

TEST_CASE("critical_point_tracking_2d_woven") {
  diy::mpi::communicator world;

  json ji = {
    {"type", "synthetic"},
    {"name", "woven"}
  };

  ftk::ndarray_stream<> stream;
  stream.configure(ji);

  ftk::critical_point_tracker_wrapper consumer;
  consumer.consume(stream);

  if (world.rank() == 0) {
    auto tracker = std::dynamic_pointer_cast<ftk::critical_point_tracker_2d_regular>( consumer.get_tracker() );
    auto trajs = tracker->get_traced_critical_points();
    auto points = tracker->get_discrete_critical_points();

    REQUIRE(trajs.size() == 48);
    // REQUIRE(points.size() == 4194); // TODO: check out if this number varies over different platform
  }
}

int main(int argc, char **argv)
{
  diy::mpi::environment env;

  Catch::Session session;
  return session.run(argc, argv);
}
