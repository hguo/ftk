#include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>

const int DW = 256, DH = 256, DT = 10;

int main(int argc, char **argv)
{
  diy::mpi::environment env;

  auto scalar = ftk::synthetic_woven_2Dt<double>(DW, DH, DT);

  ftk::critical_point_tracker_2d_regular tracker;
  tracker.set_input_scalar_field(scalar);
  // tracker.set_type_filter(ftk::CRITICAL_POINT_2D_MAXIMUM);
  tracker.update();
  tracker.write_traced_critical_points_vtk("out.vtp");

  return 0;
}
