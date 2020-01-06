#include <ftk/filters/critical_point_tracker_3d_regular.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/ndarray/conv.hh>

const int DW = 128, DH = 128, DD = 128;
const int DT = 6; // number of timesteps

int main(int argc, char **argv)
{
  diy::mpi::environment env(argc, argv);

  ftk::critical_point_tracker_3d_regular tracker(argc, argv);
  // tracker.set_domain(ftk::lattice({32, 32, 32}, {64, 64, 64}));
  // tracker.set_domain(ftk::lattice({2, 2, 2}, {DW-3, DH-3, DD-3}));
  tracker.set_domain(ftk::lattice({4, 4, 4}, {DW-5, DH-5, DD-5})); // the lattice is based on gaussian & gradient kernels
  tracker.set_array_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD}));
  tracker.set_input_array_partial(false);
  tracker.set_scalar_field_source(ftk::SOURCE_GIVEN);
  tracker.set_vector_field_source(ftk::SOURCE_DERIVED);
  tracker.set_jacobian_field_source(ftk::SOURCE_DERIVED);
  tracker.initialize();
 
  for (size_t t = 0; t < DT; t ++) {
    size_t starts[4] = {t, 0, 0, 0}, 
           sizes[4]  = {1, size_t(DD), size_t(DH), size_t(DW)};
    ftk::ndarray<double> scalar; 
    scalar.from_netcdf(argv[1], "vort", starts, sizes);
    scalar.reshape(DW, DH, DD); // reshape the 4D array to 3D

    // preconditioning
    scalar = ftk::conv3D_gaussian(scalar, 8.0/*sigma*/,
        5/*ksizex*/, 5/*ksizey*/, 5/*ksizez*/, 2/*padding*/);
    
    tracker.push_input_scalar_field(scalar);
    tracker.advance_timestep();
  }

  tracker.finalize();
  tracker.write_traced_critical_points_vtk("out.vtp");

  return 0;
}
