#include <ftk/filters/critical_point_tracker_3d_regular.hh>
#include <ftk/filters/critical_point_tracker_3d_regular_streaming.hh>
#include <ftk/filters/critical_point_tracker_3d_regular_distributed.hh>
#include <ftk/filters/critical_point_tracker_3d_regular_distributed_streaming.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>

#if FTK_HAVE_VTK
#include <ftk/geometry/points2vtk.hh>
#endif

const int DW = 128, DH = 128, DD = 128;
const int DT = 6; // number of timesteps

const std::string pattern = "/Users/hguo/workspace/data/tornado.nc/tornado*.nc";

int main(int argc, char **argv)
{
  diy::mpi::environment env(argc, argv);

  ftk::critical_point_tracker_3d_regular_distributed_streaming tracker;
  tracker.set_lb_ub({49, 48, 48, 0}, {80, 80, 80, DT});

  const auto filenames = ftk::ndarray<double>::glob(pattern);
  for (size_t t = 0; t < DT; t ++) {
    size_t starts[4] = {0, 0, 0, 0}, 
           sizes[4]  = {1, size_t(DD), size_t(DH), size_t(DW)};
    ftk::ndarray<double> u, v, w;

    // load individual vector components as a 4D array
    u.from_netcdf(filenames[t], "U", starts, sizes);
    v.from_netcdf(filenames[t], "V", starts, sizes);
    w.from_netcdf(filenames[t], "W", starts, sizes);

    ftk::ndarray<double> V; // vector field for the current timestep
    V.reshape(3, DW, DH, DD);
    for (int k = 0; k < DD; k ++)
      for (int j = 0; j < DH; j ++)
        for (int i = 0; i < DW; i ++) {
          V(0, i, j, k) = u(i, j, k, 0);
          V(1, i, j, k) = v(i, j, k, 0);
          V(2, i, j, k) = w(i, j, k, 0);
        }
    
    tracker.push_input_vector_field(V);
    tracker.advance_timestep();
  }
  tracker.update();

#if FTK_HAVE_VTK
  auto polydata = tracker.get_results_vtk();
  ftk::write_vtp("asdf11.vtp", polydata);
#endif

  return 0;
}
