#include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/ndarray/conv.hh>

const int DW = 425, DH = 880, DT = 300;

int main(int argc, char **argv)
{
  if (argc < 2) return 1;
  
  FILE *fp = fopen(argv[1], "rb");
  if (!fp) return 1;
  fseek(fp, DW*DH*sizeof(float)*200, SEEK_CUR); // skip 200 timesteps

  diy::mpi::environment env;

  ftk::critical_point_tracker_2d_regular tracker(argc, argv);
  tracker.set_domain(ftk::lattice({2, 2}, {DW-4, DH-4}));
  // tracker.set_domain(ftk::lattice({4, 4}, {DW-6, DH-6}));
  // tracker.set_domain(ftk::lattice({150, 150}, {100, 200}));
  tracker.set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
  tracker.set_input_array_partial(false);
  tracker.set_scalar_field_source(ftk::SOURCE_GIVEN);
  tracker.set_vector_field_source(ftk::SOURCE_DERIVED);
  tracker.set_jacobian_field_source(ftk::SOURCE_DERIVED);
  tracker.initialize();


  for (int k = 0; k < DT; k ++) {
    ftk::ndarray<float> scalar32({DW, DH});
    scalar32.from_binary_file(fp);

    ftk::ndarray<double> scalar;
    scalar.from_array(scalar32);
    // scalar = ftk::conv2D_gaussian(scalar, 1.0/*sigma*/, 5/*ksizex*/, 5/*ksizey*/, 2/*padding*/);

    tracker.push_input_scalar_field(scalar);
    tracker.advance_timestep();
  }

  tracker.finalize();
  tracker.write_traced_critical_points_vtk("out.vtp");

  return 0;
}
