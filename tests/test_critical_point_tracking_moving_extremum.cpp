#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include "constants.hh"
#include <ftk/filters/critical_point_tracker_wrapper.hh>

using nlohmann::json;

// #if FTK_HAVE_VTK
#if 0 // FTK_HAVE_VTK
TEST_CASE("critical_point_tracking_moving_extremum_2d") {
  auto result = track_cp2d(js_moving_extremum_2d_vti);
  diy::mpi::communicator world;
  if (world.rank() == 0)
    REQUIRE(std::get<0>(result) == 1);
}
#endif

TEST_CASE("critical_point_tracking_moving_extremum_2d_random_motion") {
  const int ncases = 10;
  const int width = 21, height = 21;
  const double x0[2] = {10.0, 10.0};
  diy::mpi::communicator world;
  
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> d{0, 1};

  for (int k = 0; k < ncases; k ++) {
    double dir[2] = {d(gen), d(gen)};
#if FTK_HAVE_MPI // if mpi is used, the exact parameters need to be the same across processes
    MPI_Bcast(dir, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    json js = js_moving_extremum_2d_synthetic;
    js["width"] = width;
    js["height"] = height;
    js["x0"] = {x0[0], x0[1]};
    js["dir"] = {dir[0], dir[1]};
    
    if (world.rank() == 0)
      std::cerr << js << std::endl;

    ftk::ndarray_stream<> stream;
    stream.configure(js);

    ftk::critical_point_tracker_wrapper consumer;
    consumer.consume(stream);

    if (world.rank() == 0) {
      auto tracker = std::dynamic_pointer_cast<ftk::critical_point_tracker_2d_regular>( consumer.get_tracker() );
      auto trajs = tracker->get_traced_critical_points();
     
      std::cerr << js << std::endl;
      // std::string f("out.txt");
      // tracker->write_traced_critical_points_text(f);

      REQUIRE(trajs.size() == 1);
      
      for (auto i = 0; i < trajs[0].size(); i ++) {
        const auto &p = trajs[0][i];
        double x = x0[0] + dir[0] * p[2], 
               y = x0[1] + dir[1] * p[2];
        // fprintf(stderr, "p={%f, %f, %f}, x={%f, %f}\n", p[0], p[1], p[2], x, y);

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
