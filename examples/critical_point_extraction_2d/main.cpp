#include <ftk/filters/extract_critical_points_2d_regular_serial.hh>
// #include <ftk/filters/extract_extrema_2d_regular_serial.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>

#if FTK_HAVE_VTK
#include <ftk/geometry/points2vtk.hh>
#endif

const int DW = 256, DH = 256;

int main(int argc, char **argv)
{
  auto scalar = ftk::synthetic_woven_2D<double>(DW, DH);
#if 0
  auto grad = hypermesh::gradient2D(scalar);
  auto hess = hypermesh::jacobian2D(grad);

  ftk::extract_critical_points_2d_regular_serial extractor;
  extractor.set_input_vector_field(grad);
  extractor.set_input_jacobian_field(hess);
  extractor.set_symmetric_jacobians(true);
  extractor.set_lb_ub({2, 2}, {DW-3, DH-3});
  extractor.set_type_filter(ftk::CRITICAL_POINT_2D_ATTRACTING ^ ftk::CRITICAL_POINT_2D_REPELLING);
  extractor.execute();
#else
  ftk::extract_critical_points_2d_regular_serial extractor;
  extractor.set_input_scalar_field(scalar);
  extractor.execute();
#endif

#if FTK_HAVE_VTK
  auto polydata = extractor.get_results_vtk();
  ftk::write_vtp("asdf.vtp", polydata);
#endif

  return 0;
}
