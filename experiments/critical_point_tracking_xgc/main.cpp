// #include <ftk/filters/critical_point_tracker_2d_regular.hh>
#include <ftk/hypermesh/simplex_2d_mesh.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/ndarray/conv.hh>

int main(int argc, char **argv)
{
  if (argc < 3) return 1; // mesh.h5, dpot.h5

#if 1 // load mesh & data from hdf5
  ftk::ndarray<int> triangles;
  triangles.from_h5(argv[1], "/cell_set[0]/node_connect_list");
  // std::cerr << triangles << std::endl; 

  ftk::ndarray<double> coords;
  coords.from_h5(argv[1], "/coordinates/values");
  // std::cerr << coords << std::endl;

  ftk::ndarray<double> dpot;
  dpot.from_h5(argv[2], "/dpot");
  dpot = dpot.transpose();
  // std::cerr << dpot << std::endl;
#endif

  ftk::simplex_2d_mesh<> m(coords, triangles);
  m.build_edges();
  m.build_smoothing_kernel(0.04);
  m.smooth_scalar_field(dpot);

#if 0
  ftk::ndarray<double> nitrz0, rxy, zxy;
  nitrz0.from_netcdf(argv[1], "NI_TRZ");
  rxy.from_netcdf(argv[1], "RXY");
  zxy.from_netcdf(argv[1], "ZXY");

  ftk::ndarray<double> coords({2, rxy.shape(0), rxy.shape(1)});
  for (size_t i = 0; i < rxy.nelem(); i ++) {
    coords[i*2] = rxy[i];
    coords[i*2+1] = zxy[i];
  }
 
#if GAUSSIAN_TIME
  auto nitrz = ftk::conv3D_gaussian(nitrz0, 
      4.0/*sigma*/, 5/*ksizex*/, 5/*ksizey*/, 5/*ksizez*/, 2/*padding*/);
      // 16.0/*sigma*/, 9/*ksizex*/, 9/*ksizey*/, 9/*ksizez*/, 4/*padding*/);
#else
  auto &nitrz = nitrz0;
#endif

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
  tracker.set_coordinates(coords);
  tracker.set_scalar_field_source(ftk::SOURCE_GIVEN);
  // tracker.set_vector_field_source(ftk::SOURCE_DERIVED);
  tracker.set_vector_field_source(ftk::SOURCE_GIVEN);
  tracker.set_jacobian_field_source(ftk::SOURCE_DERIVED);
  tracker.initialize();

  // const int ts = 500, nt = 200;
  const size_t ts = 0, nt = 704;
  for (size_t k = 0; k < nt; k ++) {
    const int t = ts+k;
    tracker.set_current_timestep(t);

    ftk::ndarray<double> scalar0 = nitrz0.slice({0, 0, ts+k}, {DW, DH, 1}); // unpreconditioned data
    scalar0.reshape({DW, DH});
#if GAUSSIAN_TIME
    ftk::ndarray<double> scalar = nitrz.slice({0, 0, ts+k}, {DW, DH, 1});
    scalar.reshape({DW, DH});
#else
    auto scalar = ftk::conv2D_gaussian(scalar0, 4.0/*sigma*/, 5/*ksizex*/, 5/*ksizey*/, 2/*padding*/);
    // auto scalar = ftk::conv2D_gaussian(scalar, 1000.0/*sigma*/, 31/*ksizex*/, 31/*ksizey*/, 15/*padding*/);
#endif

    // tracker.push_input_scalar_field(scalar);
    tracker.push_input_scalar_field(scalar0);
    tracker.push_input_vector_field(ftk::gradient2D(scalar));
    tracker.advance_timestep();
  }

  tracker.finalize();
  tracker.write_traced_critical_points("out.traj");
  tracker.write_traced_critical_points_vtk("out.vtp");
#endif
  return 0;
}
