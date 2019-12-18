#include <ftk/filters/critical_point_tracker_3d_regular.hh>
#include <ftk/filters/critical_point_tracker_3d_regular_streaming.hh>
#include <ftk/filters/critical_point_tracker_3d_regular_distributed.hh>
#include <ftk/filters/critical_point_tracker_3d_regular_distributed_streaming.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>

const int DW = 128, DH = 128, DD = 128;
const int DT = 6; // number of timesteps

int main(int argc, char **argv)
{
  diy::mpi::environment env(argc, argv);
  
  ftk::critical_point_tracker_3d_regular_distributed_streaming tracker;
  tracker.set_lb_ub({32, 32, 32, 0}, {64, 64, 64, DT});

  for (size_t t = 0; t < DT; t ++) {
    size_t starts[4] = {t, 0, 0, 0}, 
           sizes[4]  = {1, size_t(DD), size_t(DH), size_t(DW)};
    ftk::ndarray<double> scalar; 
    scalar.from_netcdf(argv[1], "vort", starts, sizes);
    scalar.reshape(DW, DH, DD); // reshape the 4D array to 3D
    
    tracker.push_input_scalar_field(scalar);
    tracker.advance_timestep();
  }
  tracker.update();
  tracker.write_traced_critical_points_vtk("out.vtp");

  return 0;
}
