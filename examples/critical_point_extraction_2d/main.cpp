#include <ftk/filters/extract_critical_points_2d_regular_serial.hh>
#include <ftk/geometry/points2vtk.hh>
#include <hypermesh/synthetic.hh>
#include <hypermesh/grad.hh>

const int DW = 256, DH = 256;

int main(int argc, char **argv)
{
  auto scalar = hypermesh::synthetic_woven_2D<double>(DW, DH);
  auto grad = hypermesh::gradient2D(scalar);
  // auto hess = hypermesh::gradient2D_vector(grad);

  ftk::extract_critical_points_2d_regular_serial extractor;
  extractor.set_input_data(grad);
  extractor.set_lb_ub({1, 1}, {DW-2, DH-2});
  extractor.execute();

  auto polydata = ftk::points2vtk(extractor.get_critical_point_coords(), 2);
  ftk::write_vtp("asdf.vtp", polydata);

  return 0;
}
