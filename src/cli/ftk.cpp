#include <fstream>
#include <mutex>
#include <set>
#include <cassert>
#include "ftk/external/cxxopts.hpp"
#include "ftk/ndarray/synthetic.hh"
#include "ftk/filters/critical_point_tracker_2d_regular.hh"
#include "ftk/filters/critical_point_tracker_3d_regular.hh"
#include "ftk/filters/json_interface.hh"
#include "ftk/filters/contour_tracker_2d_regular.hh"
#include "ftk/filters/contour_tracker_3d_regular.hh"
#include "ftk/filters/critical_line_tracker_3d_regular.hh"
#include "ftk/filters/sujudi_haimes_tracker_3d_regular.hh"
#include "ftk/filters/levy_degani_seginer_tracker_3d_regular.hh"
#include "ftk/filters/ridge_valley_tracker_3d_regular.hh"
#include "ftk/filters/tdgl_vortex_tracker_3d_regular.hh"
#include "ftk/filters/critical_line_tracker_3d_unstructured.hh"
#include "ftk/filters/xgc_blob_filament_tracker.hh"
#include "ftk/filters/xgc_blob_threshold_tracker.hh"
#include "ftk/filters/particle_tracer_regular.hh"
#include "ftk/filters/particle_tracer_mpas_ocean.hh"
#include "ftk/filters/threshold_tracker.hh"
#include "ftk/filters/streaming_filter.hh"
#include "ftk/filters/feature_curve_set_post_processor.hh"
#include "ftk/io/util.hh"
#include "ftk/io/xgc_stream.hh"
#include "ftk/io/mpas_stream.hh"
#include "ftk/ndarray.hh"
#include "ftk/ndarray/conv.hh"

using namespace ftk;
  
// global variables
std::string feature;
int ttype = 0; // feature type in integer
std::string input_pattern;
std::string output_pattern, output_type, output_format;
std::string mesh_filename;
std::string archived_intersections_filename, // archived_discrete_critical_points_filename,
  archived_traced_filename; // archived_traced_critical_points_filename;
std::string thread_backend, accelerator;
std::string type_filter_str;
int nthreads = std::thread::hardware_concurrency();
bool affinity = false, async = false;
bool verbose = false, timing = false, help = false;
int nblocks; 
bool enable_streaming_trajectories = false,
     enable_computing_degrees = false,
     disable_robust_detection = false;
int intercept_length = 2;

bool fixed_quantization_factor = false;
double quantization_factor = 0.0;

std::string post_processing_options;

size_t ntimesteps = 0, start_timestep = 0;
std::vector<std::string> input_filenames;

// contour/levelset specific
double threshold = 0.0;

// adios2 specific
std::string adios_config_file;
std::string adios_name = "BPReader";

// xgc specific
std::shared_ptr<xgc_stream> xgc_data_stream;
std::string xgc_path,
  xgc_mesh_filename, 
  xgc_ff_mesh_filename,
  xgc_bfield_filename,
  xgc_oneddiag_filename,
  xgc_units_filename,
  // xgc_augmented_mesh_filename,
  xgc_smoothing_kernel_filename,
  xgc_interpolant_filename,
  xgc_write_back_filename;
bool xgc_post_process = false, 
     xgc_torus = false, 
     xgc_use_smoothing_kernel = false,
     xgc_use_roi = false;
double xgc_smoothing_kernel_size = 0.03;
int xgc_nphi = 1, xgc_iphi = 1, xgc_vphi = 1;

// particle tracer
int pt_nsteps_per_interval = 128;
int pt_nsteps_per_checkpoint = 16;
double pt_delta_t = 1.0;
std::vector<int> pt_seed_strides;
std::vector<double> pt_seed_box;

// geo particle tracer
int ptgeo_nsteps_per_day = 0; // 1024;
int ptgeo_checkpoint_per_timestep = 0; // 
int ptgeo_checkpoint_months = 0;
int ptgeo_checkpoint_days = 0;

// mpas-o specific
std::shared_ptr<mpas_stream> mpas_data_stream;
bool geo_output = false;

// devices
std::string device_ids;
int device_buffer_size = 512; // in MB

// tracker and input stream
std::shared_ptr<tracker> mtracker;
std::shared_ptr<json_interface> wrapper;
std::shared_ptr<contour_tracker_regular> tracker_contour;
std::shared_ptr<tdgl_vortex_tracker_3d_regular> tracker_tdgl;
std::shared_ptr<critical_line_tracker> tracker_critical_line;
std::shared_ptr<critical_line_tracker_3d_regular> tracker_critical_line_regular;
std::shared_ptr<critical_line_tracker_3d_unstructured> tracker_critical_line_unstructured;
std::shared_ptr<particle_tracer> tracker_particle;
std::shared_ptr<particle_tracer_mpas_ocean> tracker_particle_mpas_ocean;
std::shared_ptr<threshold_tracker<>> tracker_threshold;
std::shared_ptr<ndarray_stream<>> stream;

nlohmann::json j_input, j_tracker;

// xgc-specific
#if 0
std::shared_ptr<simplicial_xgc_2d_mesh<>> mx2; // 2d mesh
std::shared_ptr<simplicial_xgc_3d_mesh<>> mx3, // 3d mesh, 
                                          mx30; // 3d mesh for write-backs
#endif

// input stream
static const std::string 
        str_auto("auto"),
        str_none("none"),
        str_zero("0"),
        str_two("2"),
        str_three("3"),
        str_float32("float32"),
        str_float64("float64"),
        str_adios2("adios2"),
        str_netcdf("nc"),
        str_hdf5("h5"),
        str_vti("vti"),
        str_vtp("vtp"),
        str_vtu("vtu"),
        str_scalar("scalar"),
        str_vector("vector"),
        str_text("text");

static const std::string
        str_ext_vti(".vti"), // vtkImageData
        str_ext_vtp(".vtp"), // vtkPolyData
        str_ext_vtu(".vtu"), // vtkUnstructuredGrid
        str_ext_ply(".ply"),
        str_ext_stl(".stl"),
        str_ext_netcdf(".nc"),
        str_ext_hdf5(".h5"),
        str_ext_adios2(".bp");

static const std::string
        str_critical_point_type_min("min"),
        str_critical_point_type_max("max"),
        str_critical_point_type_saddle("saddle");

static const std::set<std::string>
        set_valid_thread_backend({str_none, "pthread", "openmp", "tbb"}),
        set_valid_accelerator({str_none, "cuda", "sycl"}),
        set_valid_input_format({str_auto, str_float32, str_float64, str_netcdf, str_hdf5, str_vti, str_vtu, str_adios2}),
        set_valid_input_dimension({str_auto, str_two, str_three});

static void fatal(const cxxopts::Options &options, const std::string& str) {
  std::cerr << "FATAL: " << str << std::endl
            << options.help() << std::endl;
  exit(1);
};

static void initialize_critical_point_tracker(diy::mpi::communicator comm)
{
  j_tracker["output"] = output_pattern;
  j_tracker["output_type"] = output_type;
  j_tracker["output_format"] = output_format;
  j_tracker["intercept_length"] = intercept_length;

  j_tracker["nthreads"] = nthreads;
  j_tracker["enable_timing"] = timing;

  j_tracker["nblocks"] = std::max(comm.size(), nblocks);

  if (!mesh_filename.empty())
    j_tracker["mesh_filename"] = mesh_filename;

  if (accelerator != str_none)
    j_tracker["accelerator"] = accelerator;

  if (thread_backend != str_none)
    j_tracker["thread_backend"] = thread_backend;

  if (archived_intersections_filename.size() > 0)
    j_tracker["archived_discrete_critical_points_filename"] = archived_intersections_filename;
  
  if (archived_traced_filename.size() > 0)
    j_tracker["archived_traced_critical_points_filename"] = archived_traced_filename;
  
  if (enable_streaming_trajectories) 
    j_tracker["enable_streaming_trajectories"] = true;

  if (enable_computing_degrees)
    j_tracker["enable_computing_degrees"] = true;

  if (disable_robust_detection)
    j_tracker["enable_robust_detection"] = false;

  j_tracker["type_filter"] = type_filter_str;

  if (fixed_quantization_factor)
    j_tracker["fixed_quantization_factor"] = quantization_factor;

  if (xgc_mesh_filename.size() > 0) {
    nlohmann::json jx;
    jx["mesh_filename"] = xgc_mesh_filename;
    jx["smoothing_kernel_filename"] = xgc_smoothing_kernel_filename;
    jx["smoothing_kernel_size"] = xgc_smoothing_kernel_size;
    jx["post_process"] = xgc_post_process;
    jx["torus"] = xgc_torus;
    if (xgc_write_back_filename.size() > 0)
      jx["write_back_filename"] = xgc_write_back_filename;
    j_tracker["xgc"] = jx;
  }

  if (ttype == TRACKER_MPAS_O_CRITICAL_POINT || ttype == TRACKER_MPAS_O_PARTICLES)
    j_tracker["mpas"] = "";

  wrapper.reset(new json_interface);
  wrapper->configure(j_tracker);

  if (comm.rank() == 0) {
    // fprintf(stderr, "SUMMARY\n=============\n");
    std::cerr << "input=" << std::setw(2) << stream->get_json() << std::endl;
    std::cerr << "config=" << std::setw(2) << wrapper->get_json() << std::endl;
    // fprintf(stderr, "=============\n");
  }

  // assert(nd == 2 || nd == 3);
  // assert(nv == 1 || nv == 2 || nv == 3);
  // assert(DT > 0);
}

static void execute_critical_point_tracker(diy::mpi::communicator comm)
{
  wrapper->consume(*stream, comm);
 
  // if (!disable_post_processing)
  //    wrapper->post_process();
 
  if (!post_processing_options.empty()) {
    auto &trajs = wrapper->get_tracker()->get_traced_critical_points();
    feature_curve_set_post_processor_t pp(post_processing_options);
    pp.filter(trajs);
  }

  wrapper->write();
}

void initialize_contour_tracker(diy::mpi::communicator comm)
{
  const auto js = stream->get_json();
  const size_t nd = stream->n_dimensions(),
               DW = js["dimensions"][0], 
               DH = js["dimensions"].size() > 1 ? js["dimensions"][1].get<int>() : 0,
               DD = js["dimensions"].size() > 2 ? js["dimensions"][2].get<int>() : 0;
  const int nt = js["n_timesteps"];

  if (DD == 0) {
    tracker_contour.reset( new contour_tracker_2d_regular(comm) ); 
    tracker_contour->set_domain(lattice({0, 0}, {DW, DH}));
    tracker_contour->set_array_domain(lattice({0, 0}, {DW, DH}));
  } else {
    tracker_contour.reset(new contour_tracker_3d_regular(comm) );
    tracker_contour->set_domain(lattice({0, 0, 0}, {DW-2, DH-2, DD-2}));
    tracker_contour->set_array_domain(lattice({0, 0, 0}, {DW, DH, DD}));
  }
  
  tracker_contour->set_end_timestep(nt - 1);
  tracker_contour->set_number_of_threads(nthreads);
  tracker_contour->set_threshold(threshold);
  
  if (accelerator == "cuda")
    tracker_contour->use_accelerator(FTK_XL_CUDA);

}

void execute_contour_tracker(diy::mpi::communicator comm)
{
  tracker_contour->initialize();
  stream->set_callback([&](int k, const ndarray<double> &field_data) {
    tracker_contour->push_field_data_snapshot(field_data);
    
    if (k != 0) tracker_contour->advance_timestep();
    if (k == stream->n_timesteps() - 1) tracker_contour->update_timestep();
  });
  stream->start();
  stream->finish();
  tracker_contour->finalize();

  if (output_type == "intersections") {
    tracker_contour->write_intersections_vtp(output_pattern);
  } else if (output_type == "sliced") {
    tracker_contour->write_sliced_vtu(output_pattern);
  } else {
    tracker_contour->write_isovolume_vtu(output_pattern);
  }
}

void initialize_threshold_tracker(diy::mpi::communicator comm)
{
  tracker_threshold.reset(new threshold_tracker<>(comm));
}

void execute_threshold_tracker(diy::mpi::communicator comm)
{
  tracker_threshold->set_threshold( threshold );

  stream->set_callback([&](int k, ftk::ndarray<double> field_data) {
    fprintf(stderr, "current_timestep=%d\n", k);
    tracker_threshold->push_scalar_field_data_snapshot(field_data);
    tracker_threshold->advance_timestep();
  });

  stream->start();
  stream->finish();
  tracker_threshold->finalize();

  const auto &tg = tracker_threshold->get_tracking_graph();
  tg.generate_dot_file("dot"); // TODO
}

void initialize_xgc(diy::mpi::communicator comm)
{
  xgc_data_stream = xgc_stream::new_xgc_stream( xgc_path, comm );
  if (start_timestep == 0) start_timestep = 1;

  xgc_data_stream->set_start_timestep(start_timestep);
  xgc_data_stream->set_ntimesteps(ntimesteps);
  xgc_data_stream->set_smoothing_kernel_size( xgc_smoothing_kernel_size );
  xgc_data_stream->set_smoothing_kernel_filename( xgc_smoothing_kernel_filename );
  xgc_data_stream->set_vphi(xgc_vphi);
  xgc_data_stream->set_interpolant_filename( xgc_interpolant_filename );
  xgc_data_stream->initialize();

#if 0
  const auto js = stream->get_json();
  ntimesteps = js["n_timesteps"];

  // determine nphi and iphi
  const std::string filename0 = js["filenames"][0];
  const std::string varname = js["variables"][0];
  const auto ext = file_extension(filename0);
  
  // const auto data0 = ndarray<double>::from_file(filename0, varname);

  const auto array_nphi = ndarray<int>::from_file(filename0, "nphi");
  const auto array_iphi = ndarray<int>::from_file(filename0, "iphi");

  xgc_nphi = array_nphi[0];
  xgc_iphi = std::max(1, array_iphi[0]);

  if (xgc_path.length() > 0) {
    std::string postfix = ext == FILE_EXT_BP ? "bp" : "h5"; // TODO: check if format is bp
    xgc_mesh_filename = xgc_path + "/xgc.mesh." + postfix;
    if (!file_exists(xgc_mesh_filename)) {
      postfix = "bp";
      xgc_mesh_filename = xgc_path + "/xgc.mesh." + postfix;
    }
    if (!file_exists(xgc_mesh_filename)) {
      fatal("xgc mesh file not found.");
    }
    xgc_bfield_filename = xgc_path + "/xgc.bfield." + postfix;
    xgc_oneddiag_filename = xgc_path + "/xgc.oneddiag." + postfix;
    xgc_units_filename = xgc_path + "/units.m";
  }

  if (comm.rank() == 0) {
    fprintf(stderr, "SUMMARY\n=============\n");
    fprintf(stderr, "xgc_path=%s\n", xgc_path.c_str());
    fprintf(stderr, "xgc_mesh_filename=%s\n", xgc_mesh_filename.c_str());
    fprintf(stderr, "xgc_bfield_filename=%s\n", xgc_bfield_filename.c_str());
    fprintf(stderr, "xgc_oneddiag_filename=%s\n", xgc_oneddiag_filename.c_str());
    fprintf(stderr, "xgc_units_m=%s\n", xgc_units_filename.c_str());
    fprintf(stderr, "xgc_ff_mesh=%s\n", xgc_ff_mesh_filename.c_str());
    fprintf(stderr, "xgc_interpolant_filename=%s\n", xgc_interpolant_filename.c_str());
    fprintf(stderr, "xgc_use_smoothing_kernel=%d\n", xgc_use_smoothing_kernel);
    fprintf(stderr, "xgc_smoothing_kernel_filename=%s\n", xgc_smoothing_kernel_filename.c_str());
    fprintf(stderr, "xgc_use_roi=%d\n", xgc_use_roi);
    // fprintf(stderr, "xgc_augmented_mesh=%s\n", xgc_augmented_mesh_filename.c_str());
    fprintf(stderr, "nphi=%d, iphi=%d, vphi=%d\n", xgc_nphi, xgc_iphi, xgc_vphi);
    fprintf(stderr, "threshold=%f\n", threshold);
    fprintf(stderr, "ntimesteps=%zu\n", ntimesteps);
    fprintf(stderr, "write_back_filename=%s\n", xgc_write_back_filename.c_str());
    fprintf(stderr, "archived_intersections_filename=%s, exists=%d\n", archived_intersections_filename.c_str(), file_exists(archived_intersections_filename));
    fprintf(stderr, "archived_traced_filename=%s, exists=%d\n", archived_traced_filename.c_str(), file_exists(archived_traced_filename));
    fprintf(stderr, "output_format=%s\n", output_format.c_str());
    fprintf(stderr, "thread_backend=%s\n", thread_backend.c_str());
    fprintf(stderr, "nthreads=%d\n", nthreads);
    std::cerr << "input=" << js << std::endl;
    fprintf(stderr, "=============\n");
  }
  
  mx2 = simplicial_xgc_2d_mesh<>::from_xgc_mesh_file(xgc_mesh_filename, comm);
  // mx2 = simplicial_xgc_2d_mesh<>::from_xgc_mesh_adios2(comm, xgc_mesh_filename);
  // mx2->to_vtu("xgc_base_mesh.vtu");
 
  mx3.reset( new ftk::simplicial_xgc_3d_mesh<>(mx2, xgc_nphi, xgc_iphi, xgc_vphi) );
  
  mx2->initialize_point_locator();
  if (xgc_bfield_filename.length() > 0)
    mx2->read_bfield(xgc_bfield_filename);
  if (xgc_units_filename.length() > 0)
    mx2->read_units_m(xgc_units_filename);
  if (xgc_oneddiag_filename.length() > 0)
    mx2->read_oneddiag(xgc_oneddiag_filename);
  
  if (file_exists(archived_traced_filename))
    return; // skip interpolants, smoothing kernel

  if (xgc_use_smoothing_kernel) {
    if (file_exists(xgc_smoothing_kernel_filename))
      mx2->read_smoothing_kernel(xgc_smoothing_kernel_filename);
    else {
      mx2->build_smoothing_kernel(xgc_smoothing_kernel_size);
      if (!xgc_smoothing_kernel_filename.empty())
        mx2->write_smoothing_kernel(xgc_smoothing_kernel_filename);
    }
  }
  
  if (xgc_interpolant_filename.length() > 0) {
    if (file_exists( xgc_interpolant_filename ))
      mx3->read_interpolants( xgc_interpolant_filename );
    else {
      mx3->initialize_interpolants();
      mx3->write_interpolants( xgc_interpolant_filename );
    }
  } else 
    mx3->initialize_interpolants();
#endif
}

void initialize_xgc_blob_filament_tracker(diy::mpi::communicator comm)
{
  initialize_xgc(comm);

  std::shared_ptr<xgc_blob_filament_tracker> tracker;

  auto mx3 = xgc_data_stream->get_mx3();
  fprintf(stderr, "mx3@, n0=%zu, n1=%zu, n2=%zu, n3=%zu\n", 
      mx3->n(0), mx3->n(1), mx3->n(2), mx3->n(3));

  tracker.reset(new xgc_blob_filament_tracker(comm, xgc_data_stream->get_mx3()));
  tracker->set_device_ids(device_ids);
  tracker->set_device_buffer_size( device_buffer_size );
  tracker->set_use_roi(xgc_use_roi);
  tracker->set_end_timestep(ntimesteps - 1);
  tracker->use_thread_backend(thread_backend);
  tracker->use_accelerator(accelerator);
  tracker->set_number_of_threads(nthreads);

  if (output_type == "ordinal") {
    tracker->set_ordinal_only(true);
    tracker->set_ordinal_output_pattern(output_pattern);
  }

  tracker->initialize();

  if (archived_traced_filename.length() > 0 && file_exists(archived_traced_filename)) {
    tracker->read_surfaces(archived_traced_filename);
  } else {
    if (archived_intersections_filename.empty() || file_not_exists(archived_intersections_filename)) {
      xgc_data_stream->set_callback([&](int k, std::shared_ptr<ndarray_group> g) { // const ndarray<double> &data) {
#if 0
        auto data = g->get<double>("density");
        auto scalar = data.get_transpose();
#if 0 // FTK_HAVE_VTK
        if (xgc_write_back_filename.length()) {
          auto filename = series_filename(xgc_write_back_filename, k);
          auto mx30 = new ftk::simplicial_xgc_3d_mesh<>(mx2, xgc_nphi, xgc_iphi);
          mx30->scalar_to_vtu_slices_file(filename, "scalar", scalar);
          delete mx30;
          // tracker->get_m2()->scalar_to_xgc_slices_3d_vtu(filename, "scalar", scalar, nphi, iphi); // WIP
        }
#endif
   
        tracker->push_field_data_snapshot( scalar );
#endif

        tracker->push_field_data_snapshot(g);

        if (k != start_timestep) tracker->advance_timestep();
        if (k == stream->n_timesteps() - 1) tracker->update_timestep();
        // tracker->push_field_data_snapshot(scalar);
        // tracker->advance_timestep();
        // if (k == stream->n_timesteps() - 1) tracker->update_timestep();
      });

      while (xgc_data_stream->advance_timestep()) { }
    } else {
      fprintf(stderr, "reading archived intersections..\n");
      tracker->read_intersections_binary( archived_intersections_filename );
      tracker->set_current_timestep(stream->n_timesteps() - 1);
    }

    if ((!archived_intersections_filename.empty()) && file_not_exists(archived_intersections_filename)) {
      fprintf(stderr, "writing archived intersections..\n");
      tracker->write_intersections_binary(archived_intersections_filename);
    }
    tracker->finalize();

    if (archived_traced_filename.length() > 0) {
      fprintf(stderr, "writing archived traced intersections..\n");
      tracker->write_surfaces(archived_traced_filename, "binary");
    }
  }

  if (output_type == "ordinal") {
    // nothing to do
  }
  if (output_type == "intersections")
    tracker->write_intersections(output_pattern);
  else if (output_type == "sliced")
    tracker->write_sliced(output_pattern, output_format, true);
  else 
    tracker->write_surfaces(output_pattern, output_format, true);
}

void initialize_xgc_blob_threshold_tracker(diy::mpi::communicator comm)
{
  initialize_xgc(comm);

  std::shared_ptr<xgc_blob_threshold_tracker> tracker;
  
  // std::shared_ptr<ftk::point_locator_2d<>> locator(new ftk::point_locator_2d_quad<>(m2));
  // const double x[2] = {2.3, -0.4};
  // fprintf(stderr, "locator test: %d\n", locator->locate(x));

  tracker.reset(new xgc_blob_threshold_tracker(comm, xgc_data_stream->get_mx3()));
  tracker->set_threshold( threshold );
  tracker->initialize();

  stream->set_callback([&](int k, const ndarray<double> &data) {
    auto scalar = data.get_transpose();
    tracker->push_field_data_snapshot(scalar);

#if FTK_HAVE_VTK
    if (xgc_write_back_filename.length()) {
      auto filename = series_filename(xgc_write_back_filename, k);
      xgc_data_stream->get_mx3()->scalar_to_vtu_slices_file(filename, "scalar", scalar);
    }
#endif
#if 1 // only testing write-backs for now
    if (k != 0) tracker->advance_timestep();
    if (k == stream->n_timesteps() - 1) tracker->update_timestep();

    if (k > 0)
      tracker->write_sliced(k-1, output_pattern);
#endif      
    // mx->scalar_to_vtu_slices_file(filename, "scalar", scalar);
  });

  stream->start();
  stream->finish();

  tracker->finalize();
}

void initialize_particle_tracer_mpas_ocean(diy::mpi::communicator comm)
{
  mpas_data_stream.reset(new mpas_stream(input_pattern, comm));
  mpas_data_stream->initialize();
  if (ntimesteps != 0) {
    mpas_data_stream->set_ntimesteps(ntimesteps);
    fprintf(stderr, "overriding ntimesteps to %zu\n", ntimesteps);
  }
  // mpas_mesh_->surface_cells_to_vtu("mpas-surface.vtu");
  // exit(1);
  
  tracker_particle_mpas_ocean.reset(new particle_tracer_mpas_ocean(comm, mpas_data_stream->mesh()) );
  tracker_particle_mpas_ocean->set_number_of_threads(nthreads);
  if (accelerator == "cuda") {
    tracker_particle_mpas_ocean->use_accelerator(FTK_XL_CUDA);
    mpas_data_stream->set_c2v(false); // no need to do c2v interpolation on cpu
  }
 
  fprintf(stderr, "pt_nsteps_per_interval=%d, pt_nsteps_per_checkpoint=%d, pt_delta_t=%f\n", 
      pt_nsteps_per_interval, pt_nsteps_per_checkpoint, pt_delta_t);

  tracker_particle_mpas_ocean->set_ntimesteps(mpas_data_stream->ntimesteps);
  if (ptgeo_nsteps_per_day)
    tracker_particle_mpas_ocean->set_nsteps_per_day( ptgeo_nsteps_per_day );

  // checkpoint
  if (ptgeo_checkpoint_days)
    tracker_particle_mpas_ocean->set_checkpoint_days(ptgeo_checkpoint_days);
  else if (ptgeo_checkpoint_months)
    tracker_particle_mpas_ocean->set_checkpoint_months(ptgeo_checkpoint_months);

  // tracker_particle_mpas_ocean->set_nsteps_per_interval(pt_nsteps_per_interval);
  // tracker_particle_mpas_ocean->set_nsteps_per_checkpoint(pt_nsteps_per_checkpoint);
  tracker_particle_mpas_ocean->set_delta_t( pt_delta_t );

  tracker_particle_mpas_ocean->initialize();
  if (pt_seed_box.size() > 0)
    tracker_particle_mpas_ocean->initialize_particles_latlonz(
        pt_seed_strides[0], pt_seed_box[0], pt_seed_box[1],
        pt_seed_strides[1], pt_seed_box[2], pt_seed_box[3],
        pt_seed_strides[2], pt_seed_box[4], pt_seed_box[5]);
  else
    tracker_particle_mpas_ocean->initialize_particles_at_grid_points(pt_seed_strides);

#if 0
    ndarray<double> surfV;
    surfV.reshape(3, V.dim(2));
    surfV.set_multicomponents();
    for (size_t i = 0; i < V.dim(2); i ++)
      for (size_t k = 0; k < 3; k ++) 
        surfV(k, i) = V(k, 0, i);
    mpas_mesh_->surface_cells_to_vtu("mpas-surface.vtu", {&surfV});
    exit(1);
#endif

  mpas_data_stream->set_callback([&](int k, std::shared_ptr<ndarray_group> g) {
    tracker_particle_mpas_ocean->push_field_data_snapshot(g); // field_data);

    if (k != 0) tracker_particle_mpas_ocean->advance_timestep();
    if (k == ntimesteps - 1) tracker_particle_mpas_ocean->update_timestep();
  });

  while (mpas_data_stream->advance_timestep()) { }

  tracker_particle_mpas_ocean->finalize();
  
  if (geo_output)
    tracker_particle_mpas_ocean->write_geo_trajectories(output_pattern);
  else 
    tracker_particle_mpas_ocean->write_trajectories(output_pattern);
}

void initialize_particle_tracer(diy::mpi::communicator comm)
{
  const auto js = stream->get_json();
  const size_t nd = stream->n_dimensions(),
               DW = js["dimensions"][0], 
               DH = js["dimensions"].size() > 1 ? js["dimensions"][1].get<int>() : 0,
               DD = js["dimensions"].size() > 2 ? js["dimensions"][2].get<int>() : 0;
  const int nt = js["n_timesteps"];
  // const int nd = DD == 0 ? 2 : 3;

  std::shared_ptr<particle_tracer_regular> tr(new particle_tracer_regular(comm, nd));
  tracker_particle = tr;
  
  tr->set_ntimesteps(nt);
  tr->set_nsteps_per_interval(pt_nsteps_per_interval);
  tr->set_nsteps_per_checkpoint(pt_nsteps_per_checkpoint);
  tr->set_delta_t( pt_delta_t ); // TODO: for now, a day per timestep

  if (nd == 2) tr->set_domain(lattice({0, 0}, {DW-1, DH-1}));
  else tr->set_domain(lattice({0, 0, 0}, {DW-1, DH-1, DD-1}));
  tr->set_array_domain(lattice({0, 0}, {DW, DH}));
  tr->set_end_timestep(nt - 1);
  tr->set_number_of_threads(nthreads);

  if (comm.rank() == 0) {
    std::cerr << "input=" << std::setw(2) << stream->get_json() << std::endl;
  }

  tracker_particle->initialize();
  tracker_particle->initialize_particles_at_grid_points();

  stream->set_callback([&](int k, const ndarray<double>& field_data) {
    tracker_particle->push_field_data_snapshot("vector", field_data);

    if (k != 0) tracker_particle->advance_timestep();
    if (k == stream->n_timesteps() - 1) tracker_particle->update_timestep();
  });
  stream->start();
  stream->finish();

  tracker_particle->finalize();
  tracker_particle->write_trajectories(output_pattern);
}

void initialize_critical_line_tracker(diy::mpi::communicator comm)
{
  const auto js = stream->get_json();
  std::cerr << js << std::endl;

  if (js["format"] == "vtu") {
    auto m3 = simplicial_unstructured_3d_mesh<>::from_file(
        js["filenames"][0]);
    tracker_critical_line_unstructured.reset(new critical_line_tracker_3d_unstructured(comm, m3) );
    tracker_critical_line = tracker_critical_line_unstructured;
  } else { // regular
    const size_t nd = stream->n_dimensions(),
                 DW = js["dimensions"][0], 
                 DH = js["dimensions"].size() > 1 ? js["dimensions"][1].get<int>() : 0,
                 DD = js["dimensions"].size() > 2 ? js["dimensions"][2].get<int>() : 0;
    const int nt = js["n_timesteps"];

    if (feature == "sujudi_haimes")
      tracker_critical_line_regular.reset(new sujudi_haimes_tracker_3d_regular(comm) );
    else if (feature == "levy_degani_seginer")
      tracker_critical_line_regular.reset(new levy_degani_seginer_tracker_3d_regular(comm) );
    else if (feature == "ridge_valley")
      tracker_critical_line_regular.reset(new ridge_valley_tracker_3d_regular(comm) );
    else
      tracker_critical_line_regular.reset(new critical_line_tracker_3d_regular(comm) );
    
    tracker_critical_line_regular->set_domain(lattice({2, 2, 2}, {DW-2, DH-2, DD-2}));
    tracker_critical_line_regular->set_array_domain(lattice({0, 0, 0}, {DW, DH, DD}));
    tracker_critical_line_regular->set_end_timestep(nt - 1);
    
    if (comm.rank() == 0) {
      // fprintf(stderr, "SUMMARY\n=============\n");
      std::cerr << "input=" << std::setw(2) << stream->get_json() << std::endl;
      // fprintf(stderr, "=============\n");
    }
    
    if (accelerator == "cuda")
      tracker_critical_line_regular->use_accelerator(FTK_XL_CUDA);

    tracker_critical_line = tracker_critical_line_regular;
  }
    
  tracker_critical_line->set_number_of_threads(nthreads);
}

void execute_critical_line_tracker(diy::mpi::communicator comm)
{
  tracker_critical_line->initialize();
  stream->set_callback([&](int k, const ndarray<double> &field_data) {
    tracker_critical_line->push_field_data_snapshot(field_data);
    
    if (k != 0) tracker_critical_line->advance_timestep();
    if (k == stream->n_timesteps() - 1) tracker_critical_line->update_timestep();
  });
  stream->start();
  stream->finish();
  
  if (output_type == "intersections")
    tracker_critical_line->write_intersections(output_pattern);
  
  tracker_critical_line->finalize();
  if (output_type == "sliced")
    tracker_critical_line->write_sliced(output_pattern);
  else if (output_type == "traced")
    tracker_critical_line->write_surfaces(output_pattern, output_format);
}

void initialize_tdgl_tracker(diy::mpi::communicator comm)
{
  input_filenames = glob(input_pattern);
  if (input_filenames.empty()) 
    fatal("unable to find matching filenames.");

  if (ntimesteps != 0) input_filenames.resize(ntimesteps);

  tdgl_reader meta_reader(input_filenames[0], false); // read metadata without read actuall data
  meta_reader.read();
  const auto &meta = meta_reader.get_meta();
  
  const size_t DW = meta.dims[0], DH = meta.dims[1], DD = meta.dims[2];
  ntimesteps = input_filenames.size();
  
  if (comm.rank() == 0) {
    fprintf(stderr, "SUMMARY\n=============\n");
    fprintf(stderr, "input_pattern=%s\n", input_pattern.c_str());
    fprintf(stderr, "dims=%zu, %zu, %zu\n", DW, DH, DD);
    fprintf(stderr, "nt=%zu\n", ntimesteps);
    fprintf(stderr, "output_format=%s\n", output_format.c_str());
    fprintf(stderr, "thread_backend=%s\n", thread_backend.c_str());
    fprintf(stderr, "nthreads=%d\n", nthreads);
    fprintf(stderr, "=============\n");
  }

  tracker_tdgl.reset(new tdgl_vortex_tracker_3d_regular(comm));
  tracker_tdgl->set_domain(lattice({0, 0, 0}, {DW-2, DH-2, DD-2}));
  tracker_tdgl->set_array_domain(lattice({0, 0, 0}, {DW, DH, DD}));
  tracker_tdgl->set_end_timestep(ntimesteps - 1);
  tracker_tdgl->set_number_of_threads(nthreads);
  if (accelerator == "cuda")
    tracker_tdgl->use_accelerator(ftk::FTK_XL_CUDA);
}

void execute_tdgl_tracker(diy::mpi::communicator comm)
{
  tracker_tdgl->initialize();

  if (!archived_intersections_filename.empty()) {
    // TODO
  } else if (!archived_traced_filename.empty()) {
    tracker_tdgl->read_surfaces(archived_traced_filename);
  } else {
    for (int k = 0; k < input_filenames.size(); k ++) {
      tdgl_reader reader(input_filenames[k]);
      reader.read();

      tracker_tdgl->push_field_data_snapshot(reader.meta, 
          reader.rho, reader.phi, 
          reader.re, reader.im);
      
      if (k != 0) tracker_tdgl->advance_timestep();
      if (k == ntimesteps - 1) tracker_tdgl->update_timestep();
    }
    tracker_tdgl->finalize();
  }

  if (output_type == "intersections")
    tracker_tdgl->write_intersections(output_pattern);
  else if (output_type == "sliced")
    tracker_tdgl->write_sliced(output_pattern);
  else 
    tracker_tdgl->write_surfaces(output_pattern, output_format);
}

/////////////////
static inline nlohmann::json args_to_input_stream_json(cxxopts::ParseResult& results)
{
  using nlohmann::json;
  json j;

  if (results.count("input")) {
    const std::string input = results["input"].as<std::string>();
    if (ends_with(input, ".json")) {
      std::ifstream t(input);
      std::string str((std::istreambuf_iterator<char>(t)),
                       std::istreambuf_iterator<char>());
      t.close();
      return json::parse(str);
    }
    else 
      j["filenames"] = results["input"].as<std::string>();
  }

  if (results.count("async"))
    j["async"] = true;

  if (results.count("synthetic")) {
    j["type"] = "synthetic";
    j["name"] = results["synthetic"].as<std::string>();
  } else 
    j["type"] = "file";

  if (results.count("input-format")) j["format"] = results["input-format"].as<std::string>();
  if (results.count("dim")) j["nd"] = results["dim"].as<std::string>();
  if (results.count("bounds")) j["bounds"] = results["bounds"].as<std::string>();

  std::vector<size_t> dims;
  if (results.count("depth")) {
    dims.resize(3);
    dims[2] = results["depth"].as<size_t>();
    dims[1] = results["height"].as<size_t>();
    dims[0] = results["width"].as<size_t>();
  } else if (results.count("height")) {
    dims.resize(2);
    dims[1] = results["height"].as<size_t>();
    dims[0] = results["width"].as<size_t>();
  }
  if (dims.size())
    j["dimensions"] = dims;

  if (results.count("timesteps")) j["n_timesteps"] = results["timesteps"].as<size_t>();
  if (results.count("var")) {
    const auto var = results["var"].as<std::string>();
    const auto vars = split(var, ",");
    // if (vars.size() == 1) j["variable"] = var;
    // else if (var.size() > 1) j["variable"] = vars; 
    j["variables"] = vars;
  }

  if (results.count("temporal-smoothing-kernel")) j["temporal-smoothing-kernel"] = results["temporal-smoothing-kernel"].as<double>();
  if (results.count("temporal-smoothing-kernel-size")) j["temporal-smoothing-kernel-size"] = results["temporal-smoothing-kernel-size"].as<size_t>();
  if (results.count("spatial-smoothing-kernel")) j["spatial-smoothing-kernel"] = results["spatial-smoothing-kernel"].as<double>();
  if (results.count("spatial-smoothing-kernel-size")) j["spatial-smoothing-kernel-size"] = results["spatial-smoothing-kernel-size"].as<size_t>();
  if (results.count("perturbation")) j["perturbation"] = results["perturbation"].as<double>();

  return j;
}

///////////////////////////////
int parse_arguments(int argc, char **argv, diy::mpi::communicator comm)
{
  const int argc0 = argc;

  cxxopts::Options options(argv[0]);
  options.add_options()
    ("f,feature", "Feature type: critical_point, isosurface, tdgl_vortex, or connected_component)",
     cxxopts::value<std::string>(feature))
    ("i,input", "Input file name pattern: a single file or a series of file, e.g. 'scalar.raw', 'cm1out_000*.nc'",
     cxxopts::value<std::string>(input_pattern))
    ("input-format", "Input file format (auto|float32|float64|nc|h5|vti)", cxxopts::value<std::string>())
    ("synthetic", "Use a synthetic case (woven|double_gyre|merger) as inputs", cxxopts::value<std::string>())
    ("w,width", "Width", cxxopts::value<size_t>())
    ("h,height", "Height", cxxopts::value<size_t>())
    ("d,depth", "Depth (valid only for 3D regular grid data)", cxxopts::value<size_t>())
    ("bounds", "Bounds for resampling to regular grid data", cxxopts::value<std::string>())
    ("n,timesteps", "Number of timesteps", cxxopts::value<size_t>(ntimesteps))
    ("start-timestep", "Start timestep", cxxopts::value<size_t>(start_timestep))
    ("var", "Variable name(s), e.g. `scalar', `u,v,w'.  Valid only for NetCDF, HDF5, and VTK.", cxxopts::value<std::string>())
    ("q,quantization", "Quantization factor for SoS", cxxopts::value<double>())
    ("adios-config", "ADIOS2 config file", cxxopts::value<std::string>(adios_config_file))
    ("adios-name", "ADIOS2 I/O name", cxxopts::value<std::string>(adios_name))
    ("temporal-smoothing-kernel", "Temporal smoothing kernel bandwidth", cxxopts::value<double>())
    ("temporal-smoothing-kernel-size", "Temporal smoothing kernel size", cxxopts::value<size_t>())
    ("spatial-smoothing-kernel", "Spatial smoothing kernel bandwidth", cxxopts::value<double>())
    ("spatial-smoothing-kernel-size", "Spatial smoothing kernel size", cxxopts::value<size_t>())
    ("perturbation", "Gaussian perturbation sigma", cxxopts::value<double>())
    ("threshold", "Threshold", cxxopts::value<double>(threshold))
    ("m,mesh", "Input mesh file (will shadow arguments including width, height, depth)", cxxopts::value<std::string>(mesh_filename))
    ("nblocks", "Number of total blocks", cxxopts::value<int>(nblocks))
    // ("archived-discrete-critical-points", "Archived discrete critical points", cxxopts::value<std::string>(archived_discrete_critical_points_filename))
    ("archived-intersections", "Archived discrete intersections", cxxopts::value<std::string>(archived_intersections_filename))
    ("archived-traced", "Archived traced results", cxxopts::value<std::string>(archived_traced_filename))
    ("xgc-path", "XGC data path; will automatically read mesh, bfield, and units.m files", cxxopts::value<std::string>(xgc_path))
    ("xgc-mesh", "XGC mesh file", cxxopts::value<std::string>(xgc_mesh_filename))
    ("xgc-bfield", "XGC bfield file", cxxopts::value<std::string>(xgc_bfield_filename))
    ("xgc-ff-mesh", "XGC field following mesh file", cxxopts::value<std::string>(xgc_ff_mesh_filename))
    ("xgc-vphi", "XGC number of virtual poloidal planes", cxxopts::value<int>(xgc_vphi)->default_value("1"))
    ("xgc-roi", "XGC: extract features in ROI", cxxopts::value<bool>(xgc_use_roi))
    ("xgc-smoothing-kernel-file", "XGC: smoothing kernel file", cxxopts::value<std::string>(xgc_smoothing_kernel_filename))
    ("xgc-smoothing-kernel-size", "XGC: smoothing kernel size", cxxopts::value<double>(xgc_smoothing_kernel_size))
    ("xgc-interpolant-file", "XGC: interpolant file", cxxopts::value<std::string>(xgc_interpolant_filename))
    // ("xgc-torus", "XGC: track over poloidal planes (deprecated)", cxxopts::value<bool>(xgc_torus))
    ("xgc-write-back", "XGC: write original data back into vtu files", cxxopts::value<std::string>(xgc_write_back_filename))
    ("xgc-post-process", "XGC: enable post-processing", cxxopts::value<bool>(xgc_post_process))
    ("pt-nsteps-per-interval", "Particle tracing: Number of steps per interval", cxxopts::value<int>(pt_nsteps_per_interval))
    ("pt-nsteps-per-checkpoint", "Particle tracing: Number of steps per checkpoint", cxxopts::value<int>(pt_nsteps_per_checkpoint))
    ("pt-delta-t", "Particle tracing: Delta t per timestep", cxxopts::value<double>(pt_delta_t))
    ("pt-seed-strides", "Particle tracing: Seed strides, e.g., '4', '1,4', or '1,4,4'", cxxopts::value<std::string>())
    ("ptgeo-seeds", "Geo particle tracing: Seed in a geobox: nlat,lat0,lat1,nlon,lon0,lon1,nz,z0,z1", cxxopts::value<std::string>())
    ("ptgeo-nsteps-per-day", "Geo particle tracing: Number of steps per earth day", cxxopts::value<int>(ptgeo_nsteps_per_day))
    ("ptgeo-checkpoint-per-timestep", "Geo particle tracing: Number of checkpoints per timestep", cxxopts::value<int>(ptgeo_checkpoint_per_timestep))
    ("ptgeo-checkpoint-days", "Geo particle tracing: Number of days per checkpoint", cxxopts::value<int>(ptgeo_checkpoint_days))
    ("ptgeo-checkpoint-months", "Geo particle tracing: Number of months per checkpoint", cxxopts::value<int>(ptgeo_checkpoint_months))
    // ("xgc-augmented-mesh", "XGC: read/write augmented mesh", cxxopts::value<std::string>(xgc_augmented_mesh_filename))
    // ("mpas-data-path", "MPAS: data path", cxxopts::value<std::string>(mpas_data_path))
    ("o,output", "Output file, either one single file (e.g. out.vtp) or a pattern (e.g. out-%05d.vtp)", 
     cxxopts::value<std::string>(output_pattern))
    ("output-type", "Output type {discrete|traced|sliced|intercepted}, by default traced", 
     cxxopts::value<std::string>(output_type)->default_value("traced"))
    ("output-format", "Output format {text|vtp|vtu|ply}.  The default behavior is to automatically determine format by filename", 
     cxxopts::value<std::string>(output_format)->default_value(str_auto))
    ("geo", "Output with lonlatz coordinates.", 
     cxxopts::value<bool>(geo_output))
    ("intercept-length", "Length of intercepted outputs", 
     cxxopts::value<int>(intercept_length)->default_value("2"))
    ("type-filter", "Type filter: ane single or a combination of critical point types, e.g. `min', `max', `saddle', `min|max'",
     cxxopts::value<std::string>(type_filter_str))
    ("thread-backend", "Thread backends {pthread|openmp|tbb}",
     cxxopts::value<std::string>(thread_backend)->default_value(str_none))
    ("affinity", "Enable thread affinity", 
     cxxopts::value<bool>(affinity))
    ("nthreads", "Number of threads", 
     cxxopts::value<int>(nthreads))
    ("timing", "Enable timing", 
     cxxopts::value<bool>(timing))
    ("a,accelerator", "Accelerator {none|cuda|sycl} (experimental)",
     cxxopts::value<std::string>(accelerator)->default_value(str_none))
    ("device", "Device ID(s)", 
     cxxopts::value<std::string>(device_ids)->default_value("0"))
    ("device-buffer", "Device buffer size in MB", 
     cxxopts::value<int>(device_buffer_size)->default_value("512"))
    ("stream",  "Stream trajectories (experimental)",
     cxxopts::value<bool>(enable_streaming_trajectories))
    ("compute-degrees", "Compute degrees instead of types", 
     cxxopts::value<bool>(enable_computing_degrees)->default_value("false"))
    ("post-process", "Post process based on given options",
     cxxopts::value<std::string>(post_processing_options)->default_value(""))
    ("no-robust-detection", "Disable robust detection (faster than robust detection)",
     cxxopts::value<bool>(disable_robust_detection))
    ("async", "Asynchronous I/O", 
     cxxopts::value<bool>(async))
    ("v,verbose", "Verbose outputs", cxxopts::value<bool>(verbose))
    ("help", "Print usage", cxxopts::value<bool>(help));
  auto results = options.parse(argc, argv);

  // sanity check of arguments
  if (set_valid_accelerator.find(accelerator) == set_valid_accelerator.end())
    fatal(options, "invalid '--accelerator'");

  if (output_pattern.empty())
    fatal(options, "Missing '--output'.");

  // adios2
  if (results.count("adios-config")) { // use adios2 config file
    stream.reset(new ndarray_stream<>(adios_config_file, adios_name, comm));
  } else { // no adios2
    stream.reset(new ndarray_stream<>(comm));
  }
  
  ttype = tracker::str2tracker(feature);
  if (ttype == TRACKER_TDGL_VORTEX) { 
    // TDGL uses a different reader for now
  } else if (ttype == TRACKER_XGC_BLOB_FILAMENT || ttype == TRACKER_XGC_BLOB_THRESHOLD) {
    // xgc using a different reader for xgc filaments
  } else if (ttype == TRACKER_MPAS_O_PARTICLES) {
    // use mpas_data_stream
  } else {
    j_input = args_to_input_stream_json(results);
    stream->set_input_source_json(j_input);
  }
  if (results.count("xgc-smoothing-kernel-size") || results.count("xgc-smoothing-kernel-file"))
    xgc_use_smoothing_kernel = true;

  if (results.count("pt-seed-strides")) {
    const auto str = results["pt-seed-strides"].as<std::string>();
    const auto strs = split(str, ",");
    for (const auto s : strs)
      pt_seed_strides.push_back( std::atoi(s.c_str()) );
  }

  if (results.count("ptgeo-seeds")) {
    const auto str = results["ptgeo-seeds"].as<std::string>();
    const auto strs = split(str, ",");
    if (strs.size() < 9)
      fatal("Insufficient number for seeding geospatially.");

    pt_seed_strides.push_back( std::atoi(strs[0].c_str()) ); // nlat
    pt_seed_box.push_back( std::atof(strs[1].c_str()) );
    pt_seed_box.push_back( std::atof(strs[2].c_str()) );
    pt_seed_strides.push_back( std::atoi(strs[3].c_str()) ); // nlon
    pt_seed_box.push_back( std::atof(strs[4].c_str()) );
    pt_seed_box.push_back( std::atof(strs[5].c_str()) );
    pt_seed_strides.push_back( std::atoi(strs[6].c_str()) ); // nz
    pt_seed_box.push_back( std::atof(strs[7].c_str()) );
    pt_seed_box.push_back( std::atof(strs[8].c_str()) );
  }

  if (results.count("quantization")) {
    fixed_quantization_factor = true;
    quantization_factor = results["quantization"].as<double>();
  }

  if (ttype == TRACKER_CRITICAL_POINT)
    initialize_critical_point_tracker(comm);
  else if (ttype == TRACKER_MPAS_O_PARTICLES)
    initialize_particle_tracer_mpas_ocean(comm);
  else if (ttype == TRACKER_CRITICAL_LINE)
    initialize_critical_line_tracker(comm);
  else if (ttype == TRACKER_CONTOUR)
    initialize_contour_tracker(comm); 
  else if (ttype == TRACKER_TDGL_VORTEX)
    initialize_tdgl_tracker(comm);
  else if (ttype == TRACKER_PARTICLE)
    initialize_particle_tracer(comm);
  else if (ttype == TRACKER_XGC_BLOB_FILAMENT)
    initialize_xgc_blob_filament_tracker(comm);
  else if (ttype == TRACKER_XGC_BLOB_THRESHOLD)
    initialize_xgc_blob_threshold_tracker(comm);
  else 
    fatal(options, "missing or invalid '--feature'");

  return 0;
}

int main(int argc, char **argv)
{
  diy::mpi::environment env(argc, argv);
  diy::mpi::communicator comm; // world
  // diy::mpi::communicator comm = wcomm.split(wcomm, 20 /* a magic number as color */);
 
  parse_arguments(argc, argv, comm);

  if (ttype == TRACKER_CRITICAL_POINT || ttype == TRACKER_MPAS_O_CRITICAL_POINT)
    execute_critical_point_tracker(comm);
  else if (ttype == TRACKER_CRITICAL_LINE)
    execute_critical_line_tracker(comm);
  else if (ttype == TRACKER_CONTOUR)
    execute_contour_tracker(comm);
  else if (ttype == TRACKER_TDGL_VORTEX)
    execute_tdgl_tracker(comm);

  return 0;
}
