#include <ftk/mesh/simplicial_xgc_2d_mesh.hh>
#include <ftk/mesh/simplicial_xgc_3d_mesh.hh>
#include <ftk/io/xgc_stream.hh>
#include <ftk/features/feature_point.hh>
#include <ftk/features/feature_point_set.hh>
#include <ftk/features/feature_curve_set.hh>
#include <ftk/geometry/write_polydata.hh>
#include <ftk/geometry/points2vtk.hh>
#include <ftk/ndarray/writer.hh>
#include <ftk/filters/xgc_blob_filament_tracker.cuh>
#include <ftk/external/cxxopts.hpp>
#include <mutex>
#include <limits>

diy::mpi::environment env;

using namespace ftk;
double search_psin(
    std::shared_ptr<simplicial_xgc_2d_mesh<>> m2,
    const double psin,
    const double rz0[2], 
    const double dir[2],
    // double rz[2],
    double a, 
    double b)
{
  const double tol = 1e-8;
  while (1) {
    const double c = (a + b) * 0.5;
    const double xc[2] = {rz0[0] + dir[0] * c, rz0[1] + dir[1] * c},
                 xa[2] = {rz0[0] + dir[0] * a, rz0[1] + dir[1] * a};

    const double val = m2->eval_psin(xc) - psin;
    
    // fprintf(stderr, "a=%f, b=%f, val=%f\n", a, b, val);
    if (isnan(val)) {
      b = c;
      continue;
    } else if (val == 0.0 || (b - a) < tol) {
      // rz[0] = xc[0];
      // rz[1] = xc[1];
      return c;
    }

    if (ftk::sign(val) == ftk::sign(m2->eval_psin(xa) - psin))
      a = c;
    else 
      b = c;
  }
}

std::vector<double> initialize_seeds(
    std::shared_ptr<simplicial_xgc_2d_mesh<>> mx2,
    int nseeds_per_rake, int nrakes, 
    double psin0, double psin1, bool xpoint = false)
{
  std::vector<double> seeds(2 * nseeds_per_rake * nrakes);

  double rz0[2], rz[2];
  mx2->get_coords(0, rz0);
  ndarray<double> coords = mx2->get_coords();

  auto units = mx2->get_units();
  const double xp[2] = {units.eq_x_r, units.eq_x_z};
  const double angle0 = atan2(xp[1] - rz0[1], xp[0] - rz0[0]); //  + M_PI / 2; //  - M_PI * 3 / 2;
  
  fprintf(stderr, "rz0=%f, %f, xp=%f, %f, angle0=%f\n", rz0[0], rz0[1], xp[0], xp[1], angle0);

  const double delta_theta = 2 * M_PI / nrakes;
  for (int rake = 0; rake < nrakes; rake ++) {
    double dir[2], d0, d1;
    
    if (xpoint && rake == 0) {
      dir[0] = cos(delta_theta * rake + angle0);
      dir[1] = sin(delta_theta * rake + angle0);
      d1 = (xp[1] - rz0[1]) / dir[1];
      d0 = search_psin(mx2, psin0, rz0, dir, 0.0, d1);
    } else {
      dir[0] = cos(delta_theta * rake);
      dir[1] = sin(delta_theta * rake);
      d0 = search_psin(mx2, psin0, rz0, dir, 0.0, 2.0);
      d1 = search_psin(mx2, psin1, rz0, dir, 0.0, 2.0);
    }

    fprintf(stderr, "psin0=%f, psin1=%f, d0=%f, d1=%f\n", 
        psin0, psin1, d0, d1);

    const double offset = d0, 
                 rate = (d1 - d0) / nseeds_per_rake;

    for (int i = 0; i < nseeds_per_rake; i ++) {
      seeds[2*(rake*nseeds_per_rake+i)] = rz0[0] + dir[0] * (offset + i * rate);
      seeds[2*(rake*nseeds_per_rake+i)+1] = rz0[1] + dir[1] * (offset + i * rate);
    }
  }
  return seeds;
}

void write_apars(
    std::shared_ptr<simplicial_xgc_2d_mesh<>> mx2,
    const int np,
    const std::string& filenames, 
    const double *apars)
{
  for (int i = 0; i < np; i ++) {
    fprintf(stderr, "writing apars for plane %d.., %p\n", i, apars);
    ndarray<double> arr(apars + i * mx2->n(0), {mx2->n(0)});
    // arr.from_array(apars + i * mx2->n(0), {mx2->n(0)});
    
    auto minmax = arr.min_max();
    // fprintf(stderr, "arr_minmax=%f, %f\n", std::get<0>(minmax), std::get<1>(minmax));

    mx2->array_to_vtu( ftk::series_filename(filenames, i), 
        "apars", arr);
  }
}


void write_deltaB(
    std::shared_ptr<simplicial_xgc_2d_mesh<>> mx2,
    const int np,
    const std::string& filenames, 
    const double *deltaB)
{
  for (int i = 0; i < np; i ++) {
    fprintf(stderr, "writing deltaB for plane %d..\n", i);
    ndarray<double> arr(deltaB + i * mx2->n(0) * 3, {3, mx2->n(0)});
    arr.set_multicomponents();

    mx2->array_to_vtu( ftk::series_filename(filenames, i), 
        "deltaB", arr);
  }
}

void write_results(
    const std::string& filename,
    const double *results,
    const double *poincare_psin,
    const int nseeds_per_rake,
    const int nrakes, 
    const int nrevs)
{
  vtkSmartPointer<vtkPolyData> poly = vtkPolyData::New();
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> verts = vtkCellArray::New();

  vtkSmartPointer<vtkIntArray> id = vtkIntArray::New();
  id->SetNumberOfValues(nrevs * nrakes * nseeds_per_rake);
  id->SetName("id");

  vtkSmartPointer<vtkIntArray> rev = vtkIntArray::New();
  rev->SetNumberOfValues(nrevs * nrakes * nseeds_per_rake);
  rev->SetName("rev");

  vtkSmartPointer<vtkFloatArray> psin = vtkFloatArray::New();
  psin->SetNumberOfValues(nrevs * nrakes * nseeds_per_rake);
  psin->SetName("psin");

  vtkSmartPointer<vtkFloatArray> psin0 = vtkFloatArray::New();
  psin0->SetNumberOfValues(nrevs * nrakes * nseeds_per_rake);
  psin0->SetName("psin0");

  vtkSmartPointer<vtkFloatArray> delta_psin = vtkFloatArray::New();
  delta_psin->SetNumberOfValues(nrevs * nrakes * nseeds_per_rake);
  delta_psin->SetName("delta_psin");
  
  vtkSmartPointer<vtkFloatArray> offset_psin0 = vtkFloatArray::New();
  offset_psin0->SetNumberOfValues(nrevs * nrakes * nseeds_per_rake);
  offset_psin0->SetName("offset_psin0");

  std::vector<int> arr_rev;

  for (int i = 0; i < nrevs; i ++) {
    for (int j = 0; j < nrakes; j ++) {
      for (int k = 0; k < nseeds_per_rake; k ++) {
        const int idx = i * nrakes * nseeds_per_rake + j * nseeds_per_rake + k;
        const int idx0 = j* nseeds_per_rake + k;
        const int idx_last = i == 0 ? idx0 : 
          ((i-1) * nrakes * nseeds_per_rake + j * nseeds_per_rake + k);
        // fprintf(stderr, "idx=%d\n", idx);

        double p[3] = {results[idx*2], results[idx*2+1], 0.0};
        double psin_ = poincare_psin[idx];
        double psin0_ = poincare_psin[idx0];
        double psin_last = poincare_psin[idx_last];
        double delta_psin_ = psin_last - psin_;
        double offset_psin0_ = psin_ - psin0_;
        if (p[0] > 1e37) {
          p[0] = p[1] = p[2] = std::nan("");
          psin_ = std::nan("");
          psin0_ = std::nan("");
          delta_psin_ = std::nan("");
          offset_psin0_ = std::nan("");
        }

        vtkIdType pid = points->InsertNextPoint(p);
        verts->InsertNextCell(1, &pid);
        
        rev->SetValue(idx, i);
        psin->SetValue(idx, psin_);
        psin0->SetValue(idx, psin0_);
        delta_psin->SetValue(idx, delta_psin_);
        offset_psin0->SetValue(idx, offset_psin0_);
        id->SetValue(idx, idx0);
      }
    }
  }

  poly->SetPoints(points);
  poly->SetVerts(verts);
  poly->GetPointData()->AddArray(id);
  poly->GetPointData()->AddArray(psin);
  poly->GetPointData()->AddArray(psin0);
  poly->GetPointData()->AddArray(delta_psin);
  poly->GetPointData()->AddArray(offset_psin0);
  poly->GetPointData()->AddArray(rev);

  write_polydata(filename, poly);
}

int main(int argc, char **argv)
{
  int device;
  std::string path, input, output, deltaB_filenames, apars_filenames;
  double psin0 = 0.8, psin1 = 1.03;
  int vphi = 64;
  int nseeds_per_rake = 500, nrakes = 6;
  int nrevs = 3000;
  bool trace_static = false;
  bool trace_reverse = false;
  bool xpoint = false;

  cxxopts::Options options(argv[0]);
  options.add_options()
    ("p,path", "XGC data path; will automatically read mesh, bfield, and units.m files", cxxopts::value<std::string>(path))
    ("i,input", "input", cxxopts::value<std::string>(input))
    ("d,device", "device id", cxxopts::value<int>(device)->default_value("0"))
    ("static", "trace with static magnetic field", cxxopts::value<bool>(trace_static))
    ("reverse", "trace in reverse direction", cxxopts::value<bool>(trace_reverse))
    ("psin0", "psin0", cxxopts::value<double>(psin0)->default_value("0.8"))
    ("psin1", "psin1", cxxopts::value<double>(psin1)->default_value("1.03"))
    ("v,vphi", "number of virtual poloidal planes", cxxopts::value<int>(vphi)->default_value("64"))
    ("n,rev", "revolutions", cxxopts::value<int>(nrevs)->default_value("3000"))
    ("s,seeds", "number of seeds per rake", cxxopts::value<int>(nseeds_per_rake)->default_value("500"))
    ("r,rakes", "number of rakes", cxxopts::value<int>(nrakes)->default_value("6"))
    ("x,xpoint", "include x-point", cxxopts::value<bool>(xpoint))
    ("o,output", "output file", cxxopts::value<std::string>(output))
    ("dump-deltaB", "dump deltaB for debugging", cxxopts::value<std::string>(deltaB_filenames))
    ("dump-apars", "dump apars for debugging", cxxopts::value<std::string>(apars_filenames));
  auto parse_results = options.parse(argc, argv);

  fprintf(stderr, "nseeds_per_rake=%d, nrakes=%d, nrevs=%d\n", 
      nseeds_per_rake, nrakes, nrevs);

  auto xs = xgc_stream::new_xgc_stream(path);

  xs->set_enable_initialize_smoothing_kernel( false );
  xs->set_enable_initialize_interpolants( false );
  xs->set_vphi( vphi );
  xs->initialize();

  auto mx2 = xs->get_m2(); // simplicial_xgc_2d_mesh<>::from_xgc_mesh_file(argv[1]);
  auto mx3 = xs->get_mx3();

  // const auto psifield = mx2->get_psifield();
  // ndarray<double> psinfield = psifield * (1.0 / mx2->get_units().psi_x);
  // auto psin_minmax = psinfield.min_max();
  // fprintf(stderr, "psin_min=%f, psin_max=%f\n", std::get<0>(psin_minmax), std::get<1>(psin_minmax));

  const int nseeds = nseeds_per_rake * nrakes;
  auto seeds = initialize_seeds(
      mx2, nseeds_per_rake, nrakes, psin0, psin1, xpoint);

  // mx2->initialize_point_locator();
  // mx2->read_bfield(argv[2]);
  mx2->derive_curl_bfield0();
 
#if 0
  const int nphi = apars.dim(1), iphi = 1;
  fprintf(stderr, "nphi=%d, iphi=%d, vphi=%d\n", nphi, iphi, vphi);
  // const int nphi = 16, iphi = 1, vphi = 16; 
  auto mx3 = new simplicial_xgc_3d_mesh<>(mx2, nphi, iphi, vphi);
#endif

#if 0
  ftk::ndarray<double> apars_upsample = mx3->interpolate(apars);
  apars_upsample.set_multicomponents();
  std::cerr << apars_upsample.shape() << std::endl;
#endif

  const ndarray<double>& bfield = mx2->get_bfield();
  const ndarray<double>& bfield0 = mx2->get_bfield0();
  const ndarray<double>& curl_bfield0 = mx2->get_curl_bfield0();

  const int n0 = mx2->n(0); // apars.dim(0);
  const int nphi = mx3->get_nphi(), iphi = mx3->get_iphi();

  xft_ctx_t *ctx;
  xft_create_poincare_ctx(&ctx, nseeds, nrevs, device);
  
  xft_load_mesh(ctx, nphi, iphi, vphi, 
      mx2->n(0), mx2->n(1), mx2->n(2), 
      mx2->get_coords().data(), 
      mx2->get_edges().data(),
      mx2->get_triangles().data());

  xft_load_magnetic_field(ctx, bfield.data(), bfield0.data(), curl_bfield0.data());

  xft_load_vertex_triangles(ctx, mx2->get_vertex_triangles());

  xft_load_psin(ctx, mx2->get_psinfield().data());

  const auto bvh = std::static_pointer_cast<point_locator_2d_quad<>>(mx2->get_locator())->to_bvh();
  xft_load_bvh(ctx, bvh);

  if (!trace_static) {
#if 1
    mx3->initialize_interpolants_cached();
    xft_load_interpolants(ctx, mx3->get_interpolants());
#else
    xft_derive_interpolants(ctx);
#endif

    ftk::ndarray<double> apars;
    apars.read_bp(input, "apars");
    
    // auto minmax = apars.min_max();
    // fprintf(stderr, "apar_minmax=%f, %f\n", std::get<0>(minmax), std::get<1>(minmax));
    // std::cerr << apars.shape() << std::endl;
    // write_apars(mx2, nphi, apars_filenames, apars.data()); // diagnosis
    // exit(1);
    
    // if (apars_filenames.size() > 0)
    //   ctx->retrieve_apars_upsample = true;

    xft_load_apars(ctx, apars.data());

    if (apars_filenames.size() > 0) {
      write_apars(mx2, nphi * vphi, apars_filenames, ctx->h_apars_upsample);
      xft_destroy_ctx(&ctx);
      return 0;
    }

#if 0
    if (deltaB_filenames.size() > 0) { // dump deltaB field
      assert(false); // TODO
      xft_retrieve_deltaB(ctx);
      // write_deltaB(mx2, nphi * vphi, deltaB_filenames, ctx->h_deltaB);
      xft_destroy_ctx(&ctx);
      exit(0);
    }
#endif
  }

  xft_compute_poincare_plot(ctx, &seeds[0], 
      trace_static, 
      trace_reverse ? -1 : 1);

  xft_compute_poincare_psin(ctx);

  fprintf(stderr, "finalizing...\n");
  const double *plot = (const double*)(ctx->hcps);
#if 0
  std::vector<double> results;
  // std::vector<double> plot(plot_raw, plot_raw + nseeds*nrevs*2);
  for (int j = 0; j < nseeds; j ++) {
    for (int i = 0; i < nrevs; i ++) {
      const size_t offset = i*nseeds + j;
      const double r = plot[offset*2], z = plot[offset*2+1];
      if (r < 1e37) {
        results.push_back(r);
        results.push_back(z);
      }
      // fprintf(stderr, "%f, %f\n", plot[offset*2], plot[offset*2+1]);
    }
  }
#endif

  fprintf(stderr, "writing...\n");
  write_results(output, plot, ctx->h_poincare_psin, 
      nseeds_per_rake, nrakes, nrevs);
  
  fprintf(stderr, "exiting...\n");
  xft_destroy_ctx(&ctx);

  return 0;
}

