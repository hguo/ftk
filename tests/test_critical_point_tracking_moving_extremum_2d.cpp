#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include "constants.hh"
#include <ftk/filters/critical_point_tracker_wrapper.hh>

#if FTK_HAVE_VTK
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#endif

using nlohmann::json;

// #if FTK_HAVE_VTK
#if 0 // FTK_HAVE_VTK
TEST_CASE("critical_point_tracking_moving_extremum_2d") {
  auto result = track_cp2d(js_moving_extremum_2d_vti);
  diy::mpi::communicator comm;
  if (comm.rank() == 0)
    REQUIRE(std::get<0>(result) == 1);
}
#endif

TEST_CASE("critical_point_tracking_moving_extremum_2d_random_motion") {
  const int ncases = 10;
  const int width = 21, height = 21;
  const double x0[2] = {10.0, 10.0};
  
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> d{0, 1};

  diy::mpi::communicator comm;
  
  for (int k = 0; k < ncases; k ++) {
    int root = rand() % comm.size();
    diy::mpi::broadcast(comm, root, 0);

    std::vector<double> dir = {d(gen), d(gen)};
    diy::mpi::broadcast(comm, dir, root);

    json js = js_moving_extremum_2d_synthetic;
    js["dimensions"] = {width, height};
    js["x0"] = {x0[0], x0[1]};
    js["dir"] = {dir[0], dir[1]};
   
    if (comm.rank() == root)
      std::cerr << js << std::endl;

    ftk::ndarray_stream<> stream;
    stream.configure(js);
    
    json jc;
    jc["root_proc"] = root;
    // fprintf(stderr, "root=%d\n", root);

    ftk::critical_point_tracker_wrapper consumer;
    consumer.configure(jc);
    consumer.consume(stream);

    if (comm.rank() == root) {
      auto tracker = std::dynamic_pointer_cast<ftk::critical_point_tracker_2d_regular>( consumer.get_tracker() );
      auto trajs = tracker->get_traced_critical_points();
     
      std::cerr << js << std::endl;
      // std::string f("out.txt");
      
#if FTK_HAVE_VTK 
      // tracker->write_traced_critical_points_text(f);
      std::string filename_vtk = "moving_extremum_2d-" + std::to_string(k) + ".vtp";
      tracker->write_traced_critical_points_vtk(filename_vtk);

      vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkXMLPolyDataReader::New();
      reader->SetFileName(filename_vtk.c_str());
      reader->Update();
      vtkSmartPointer<vtkPolyData> poly = reader->GetOutput();
      REQUIRE(poly->GetNumberOfCells() == 1);
#endif

      REQUIRE(trajs.size() == 1);
      
      for (auto i = 0; i < trajs[0].size(); i ++) {
        const auto &p = trajs[0][i];
        double x = x0[0] + dir[0] * p.t, 
               y = x0[1] + dir[1] * p.t;
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
