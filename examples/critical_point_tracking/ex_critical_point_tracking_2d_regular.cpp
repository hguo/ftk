#include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/ndarray/conv.hh>

// This example tracks critical points from a synthetic data.
// If FTK is built with VTK, the output trajectories are 
// stored in out.vti, which can be opened and visualized by ParaView.

const int DW = 32, DH = 32, DT = 10;

int main(int argc, char **argv)
{
  diy::mpi::environment env;

  ftk::critical_point_tracker_2d_regular tracker;
  tracker.set_domain(ftk::lattice({2, 2}, {DW-3, DH-3}));
  // tracker.set_domain(ftk::lattice({4, 4}, {DW-6, DH-6}));
  tracker.set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
  tracker.set_input_array_partial(false);
  tracker.set_scalar_field_source(ftk::SOURCE_GIVEN);
  tracker.set_vector_field_source(ftk::SOURCE_DERIVED);
  tracker.set_jacobian_field_source(ftk::SOURCE_DERIVED);
  // tracker.set_type_filter(ftk::CRITICAL_POINT_2D_MAXIMUM);
  tracker.initialize();

  // FILE *fp = fopen("out.raw", "wb");
  for (int k = 0; k < DT; k ++) {
    auto scalar = ftk::synthetic_woven_2D<double>(DW, DH, double(k) / (DT - 1));
    // auto scalar = ftk::synthetic_merger_2D<double>(DW, DH, double(k) / (DT - 1) * 20);
    // scalar = ftk::conv2D_gaussian(scalar, 5.0/*sigma*/, 5/*ksizex*/, 5/*ksizey*/, 2/*padding*/);
    // scalar.to_binary_file(fp);

    tracker.push_scalar_field_snapshot(scalar);
    if (k != 0)
      tracker.advance_timestep();
  }
  // fclose(fp);

  // ftk::ndarray<double> scalar_all;
  // scalar_all.reshape(DW, DH, DT);
  // scalar_all.from_binary_file("out.raw");
  // scalar_all.to_scalar_vtk_image_data_file("out.vti");

  tracker.finalize();
  tracker.write_traced_critical_points_vtk("out.vtp");
  // tracker.write_discrete_critical_points_vtk("out.vtp");

  return 0;
}
