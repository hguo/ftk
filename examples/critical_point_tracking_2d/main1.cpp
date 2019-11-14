#include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>

#if FTK_HAVE_VTK
#include <ftk/geometry/points2vtk.hh>
#endif

const int DW = 256, DH = 256, DT = 10;

int main(int argc, char **argv)
{
  auto scalar = ftk::synthetic_woven_2Dt<double>(DW, DH, DT);

  ftk::critical_point_tracker_2d_regular tracker;
  tracker.set_input_scalar_field(scalar);
  // tracker.set_type_filter(ftk::CRITICAL_POINT_2D_MAXIMUM);
  tracker.update();

#if 1 // FTK_HAVE_VTK
  auto polydata = tracker.get_results_vtk();
  ftk::write_vtp("asdf1.vtp", polydata);
#endif

  return 0;
}
