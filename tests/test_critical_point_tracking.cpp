#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <ftk/filters/critical_point_tracker_wrapper.hh>

using nlohmann::json;

TEST_CASE("critical_point_tracking_2d_woven") {
  json ji = {
    {"type", "synthetic"},
    {"name", "woven"}
  };

  ftk::ndarray_stream<> stream;
  stream.configure(ji);

  ftk::critical_point_tracker_wrapper consumer;
  consumer.consume(stream);

  auto tracker = std::dynamic_pointer_cast<ftk::critical_point_tracker_2d_regular>( consumer.get_tracker() );
  auto trajs = tracker->get_traced_critical_points();
  auto points = tracker->get_discrete_critical_points();

  REQUIRE(trajs.size() == 48);
  // REQUIRE(points.size() == 4194); // TODO: this number seems to vary.  need to find out the reason
}

int main(int argc, char **argv)
{
  Catch::Session session;
  return session.run(argc, argv);
}
