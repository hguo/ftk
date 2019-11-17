#include <ftk/filters/critical_point_tracker_3d_regular.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>

#if FTK_HAVE_VTK
#include <ftk/geometry/points2vtk.hh>
#endif

const int DW = 128, DH = 128, DD = 128;// the dimensionality of the data is DW*DH
const int DT = 6; // number of timesteps

int main(int argc, char **argv)
{
  diy::mpi::environment env(argc, argv);

  size_t starts[4] = {0, 0, 0, 0}, 
         sizes[4]  = {size_t(DT), size_t(DD), size_t(DH), size_t(DW)};

  ftk::ndarray<float> scalar_f;
  scalar_f.reshape(DW, DH, DD, DT);
  scalar_f.from_netcdf(argv[1], "vort", starts, sizes);
  
  ftk::ndarray<double> scalar;
  scalar.from_array(scalar_f);

  ftk::critical_point_tracker_3d_regular tracker;
  tracker.set_input_scalar_field(scalar);
  // tracker.set_type_filter(ftk::CRITICAL_POINT_3D_MAXIMUM);
  tracker.set_lb_ub({32, 32, 32, 0}, {64, 64, 64, DT-1});
  tracker.update();

#if FTK_HAVE_VTK
  auto polydata = tracker.get_results_vtk();
  ftk::write_vtp("asdf6.vtp", polydata);
#endif

  return 0;
}
