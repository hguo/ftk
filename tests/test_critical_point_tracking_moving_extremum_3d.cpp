#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include "constants.hh"

#if FTK_HAVE_VTK
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#endif

using nlohmann::json;
using Catch::Approx;

TEST_CASE("critical_point_tracking_moving_extremum_3d_random_motion") {
  const int ncases = 3;
  const int width = 21, height = 21, depth = 21;
  const double x0[3] = {10.0, 10.0, 10.0};
  // const double x0[3] = {10.05, 10.03, 10.1};
  
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> d{0, 1};

  diy::mpi::communicator comm;
  
  for (int k = 0; k < ncases; k ++) {
    int root = rand() % comm.size();
    diy::mpi::broadcast(comm, root, 0);

    std::vector<double> dir = {d(gen), d(gen), d(gen)};
    diy::mpi::broadcast(comm, dir, root);

    json js = js_moving_extremum_3d_synthetic;
    js["dimensions"] = {width, height, depth};
    js["x0"] = {x0[0], x0[1], x0[2]};
    js["dir"] = {dir[0], dir[1], dir[2]};
    js["n_timesteps"] = 10;
   
    if (comm.rank() == root)
      std::cerr << js << std::endl;

    ftk::ndarray_stream<> stream;
    stream.configure(js);
    
    json jc;
    jc["root_proc"] = root;
    // fprintf(stderr, "root=%d\n", root);

    ftk::json_interface consumer;
    consumer.configure(jc);
    consumer.consume(stream);
    consumer.post_process();

    if (comm.rank() == root) {
      auto tracker = std::dynamic_pointer_cast<ftk::critical_point_tracker_3d_regular>( consumer.get_tracker() );
      auto trajs = tracker->get_traced_critical_points();
     
      // std::cerr << js << std::endl;
      // std::string f("out.txt");
      
#if 0 // FTK_HAVE_VTK 
      // tracker->write_traced_critical_points_text(f);
      std::string filename_vtk = "moving_extremum_3d-" + std::to_string(k) + ".vtp";
      tracker->write_traced_critical_points_vtk(filename_vtk);

      vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkXMLPolyDataReader::New();
      reader->SetFileName(filename_vtk.c_str());
      reader->Update();
      vtkSmartPointer<vtkPolyData> poly = reader->GetOutput();
      REQUIRE(poly->GetNumberOfCells() == 1);
#endif

      REQUIRE(trajs.size() == 1);
      
      for (const auto &kv : trajs) {
        const auto &traj = kv.second;
        for (auto i = 0; i < traj.size(); i ++) {
          const auto &p = traj[i];
          double x = x0[0] + dir[0] * p[3], 
                 y = x0[1] + dir[1] * p[3], 
                 z = x0[2] + dir[2] * p[3];
          fprintf(stderr, "p={%f, %f, %f, %f}, x={%f, %f, %f}\n", p[0], p[1], p[2], p[3], x, y, z);

          REQUIRE(p[0] == Approx(x).epsilon(0.01));
          REQUIRE(p[1] == Approx(y).epsilon(0.01));
          REQUIRE(p[2] == Approx(z).epsilon(0.01));
        }
      }
    }
  }
}

#include "main.hh"
