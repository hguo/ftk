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

  size_t starts[4] = {0, 0, 0, 0}, 
         sizes[4]  = {size_t(DT), size_t(DD), size_t(DH), size_t(DW)};

  ftk::ndarray<double> scalar;
  scalar.reshape(DW, DH, DD, DT);
  scalar.from_netcdf(argv[1], "vort", starts, sizes);

  ftk::critical_point_tracker_3d_regular tracker(argc, argv);
  // tracker.use_accelerator(ftk::FTK_XL_CUDA);
  tracker.set_input_scalar_field(scalar);
  // tracker.set_type_filter(ftk::CRITICAL_POINT_3D_MAXIMUM);
  tracker.set_lb_ub({32, 32, 32, 0}, {64, 64, 64, DT-1});
  tracker.update();
  tracker.write_traced_critical_points_vtk("out.vtp");

  return 0;
}
