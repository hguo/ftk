#include <ftk/filters/critical_point_tracker_2d_regular_streaming.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>

const int DW = 256, DH = 256, DT = 10;

int main(int argc, char **argv)
{
  diy::mpi::environment env;

  ftk::critical_point_tracker_2d_regular_streaming tracker;
  // tracker.set_type_filter(ftk::CRITICAL_POINT_2D_MAXIMUM);
 
  for (int k = 0; k < DT; k ++) {
    auto scalar = ftk::synthetic_woven_2D<double>(DW, DH, double(k) / (DT - 1));
    tracker.push_input_scalar_field(scalar);
    tracker.advance_timestep();
  }

  tracker.update();
  tracker.write_traced_critical_points_vtk("out.vtp");

  return 0;
}
