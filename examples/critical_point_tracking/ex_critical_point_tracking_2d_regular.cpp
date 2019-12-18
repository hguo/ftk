#include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>

const int DW = 256, DH = 256, DT = 10;

int main(int argc, char **argv)
{
  diy::mpi::environment env;

  ftk::critical_point_tracker_2d_regular tracker;
  tracker.set_domain(ftk::lattice({2, 2}, {DW-4, DH-4}));
  tracker.set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
  tracker.set_input_array_partial(false);
  tracker.set_scalar_field_source(ftk::SOURCE_GIVEN);
  tracker.set_vector_field_source(ftk::SOURCE_DERIVED);
  tracker.set_jacobian_field_source(ftk::SOURCE_DERIVED);
  tracker.initialize();
 
  for (int k = 0; k < DT; k ++) {
    fprintf(stderr, "k=%d\n", k);
    auto scalar = ftk::synthetic_woven_2D<double>(DW, DH, double(k) / (DT - 1));
    tracker.push_input_scalar_field(scalar);
    tracker.advance_timestep();
  }

  tracker.finalize();
  tracker.write_traced_critical_points_vtk("out.vtp");

  return 0;
}
