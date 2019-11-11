#include <ftk/filters/critical_point_tracker_2d_regular_distributed.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>

const int DW = 256, DH = 256, DT = 10;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  auto scalar = ftk::synthetic_woven_2Dt<double>(DW, DH, DT);
  
  ftk::critical_point_tracker_2d_regular_distributed tracker;

  tracker.set_input_scalar_field(scalar);
  tracker.update();

  MPI_Finalize();
  return 0;
}
