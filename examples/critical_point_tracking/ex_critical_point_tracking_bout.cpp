#include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/ndarray/conv.hh>

#if 1
int main(int argc, char **argv)
{
  if (argc < 2) return 1;

  ftk::ndarray<double> nitrz;
  nitrz.from_netcdf(argv[1], "NI_TRZ");
    
  nitrz  = ftk::conv3D_gaussian(nitrz, 
      8.0/*sigma*/, 5/*ksizex*/, 5/*ksizey*/, 5/*ksizez*/, 2/*padding*/);
      // 16.0/*sigma*/, 9/*ksizex*/, 9/*ksizey*/, 9/*ksizez*/, 4/*padding*/);

  const size_t DW = nitrz.shape(0), 
               DH = nitrz.shape(1),
               DT = nitrz.shape(2);

  fprintf(stderr, "DW=%lu, DH=%lu, DT=%lu\n", DW, DH, DT);

  diy::mpi::environment env;

  ftk::critical_point_tracker_2d_regular tracker(argc, argv);
  // tracker.set_domain(ftk::lattice({2, 2}, {DW-4, DH-4}));
  // tracker.set_domain(ftk::lattice({4, 4}, {DW-6, DH-6}));
  // tracker.set_domain(ftk::lattice({2300, 4}, {200, 250}));
  tracker.set_domain(ftk::lattice({2300, 150}, {200, 30}));
  tracker.set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
  tracker.set_input_array_partial(false);
  tracker.set_scalar_field_source(ftk::SOURCE_GIVEN);
  tracker.set_vector_field_source(ftk::SOURCE_DERIVED);
  tracker.set_jacobian_field_source(ftk::SOURCE_DERIVED);
  tracker.initialize();

  // const int ts = 500, nt = 200;
  const size_t ts = 0, nt = 704;
  for (size_t k = 0; k < nt; k ++) {
    const int t = ts+k;
    tracker.set_current_timestep(t);
    ftk::ndarray<double> scalar = nitrz.slice({0, 0, ts+k}, {DW, DH, 1});
    scalar.reshape({DW, DH});
    // auto scalar1 = ftk::conv2D_gaussian(scalar, 4.0/*sigma*/, 5/*ksizex*/, 5/*ksizey*/, 2/*padding*/);
    // auto scalar1 = ftk::conv2D_gaussian(scalar, 1000.0/*sigma*/, 31/*ksizex*/, 31/*ksizey*/, 15/*padding*/);

    tracker.push_input_scalar_field(scalar);
    tracker.advance_timestep();
  }

  tracker.finalize();
  tracker.write_traced_critical_points("out.traj");
  tracker.write_traced_critical_points_vtk("out.vtp");

  return 0;
}
#endif

#if 0
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
#endif
