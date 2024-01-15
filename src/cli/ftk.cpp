#include <iostream>
#include <ndarray/ndarray_group_stream.hh>
#include <ftk/external/cxxopts.hpp>
#include <ftk/filters/particle_tracer_mpas_ocean.hh>

// global variables
std::string feature;
int ttype = 0; // feature type in integer
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

// input stream
std::string input_pattern;
std::string input_prefix;
std::string &stream_yaml_filename(input_pattern); // alias
std::string mesh_filename;
size_t ntimesteps = 0, start_timestep = 0;

// devices
std::string device_ids;
int device_buffer_size = 512; // in MB

// adios2 specific
std::string adios_config_file;
std::string adios_name = "BPReader";

// output
std::string output_pattern, output_type, output_format;
std::string archived_intersections_filename, // archived_discrete_critical_points_filename,
  archived_traced_filename; // archived_traced_critical_points_filename;

// contour/levelset specific
double threshold = 0.0;

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
bool geo_output = false;

// xgc specific
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

// input stream specific
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
 
// input stream
std::shared_ptr<ftk::stream> stream;

// trackers
std::shared_ptr<ftk::particle_tracer_mpas_ocean> tracker_particle_mpas_ocean;

bool parse_arguments(int argc, char **argv)
{
  const int argc0 = argc;

  cxxopts::Options options(argv[0]);
  options.add_options()
    ("f,feature", "Feature type: critical_point, isosurface, tdgl_vortex, or connected_component)",
     cxxopts::value<std::string>(feature))
    ("i,input", "Input: a yaml file, a single file, or a series of file, e.g. 'scalar.raw', 'cm1out_000*.nc'",
     cxxopts::value<std::string>(input_pattern))
    ("p,prefix", "Input path prefix override", cxxopts::value<std::string>(input_prefix))
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
  // options.parse_positional({"input"});

  if (!results.count("input") || results.count("help")) {
    std::cerr << options.help() << std::endl;
    exit(0);
  }

  fprintf(stderr, "input stream yaml file: %s\n", stream_yaml_filename.c_str());

  stream.reset(new ftk::stream);

  if (!input_prefix.empty())
    stream->set_path_prefix( input_prefix );

  stream->parse_yaml(stream_yaml_filename);

  return 0;
}

void initialize_particle_tracer_mpas_ocean(diy::mpi::communicator comm)
{
  auto gs = stream->read_static(); // static array group
  gs->print_info(std::cerr);

  std::shared_ptr<ftk::mpas_mesh<>> mesh(new ftk::mpas_mesh<>(gs));

  tracker_particle_mpas_ocean.reset(new ftk::particle_tracer_mpas_ocean(comm, mesh) );
#if 0
  tracker_particle_mpas_ocean->set_number_of_threads(nthreads);
  if (accelerator == "cuda") {
    tracker_particle_mpas_ocean->use_accelerator(FTK_XL_CUDA);
    mpas_data_stream->set_c2v(false); // no need to do c2v interpolation on cpu
    mpas_data_stream->set_e2c(false);
  }
#endif
 
  fprintf(stderr, "pt_nsteps_per_interval=%d, pt_nsteps_per_checkpoint=%d, pt_delta_t=%f\n", 
      pt_nsteps_per_interval, pt_nsteps_per_checkpoint, pt_delta_t);

  tracker_particle_mpas_ocean->set_ntimesteps(stream->total_timesteps());
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

  for (int i = 0; i < stream->total_timesteps(); i ++)  {
    auto g = stream->read(i);
    g->print_info(std::cerr);
    
    tracker_particle_mpas_ocean->push_field_data_snapshot(g); // field_data);
    
    if (i != 0) tracker_particle_mpas_ocean->advance_timestep();
    if (i == stream->total_timesteps() - 1) tracker_particle_mpas_ocean->update_timestep();
  }

  tracker_particle_mpas_ocean->finalize();
  
  if (geo_output)
    tracker_particle_mpas_ocean->write_geo_trajectories(output_pattern);
  else 
    tracker_particle_mpas_ocean->write_trajectories(output_pattern);
}

int main(int argc, char **argv)
{
#if NDARRAY_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  parse_arguments(argc, argv);

  initialize_particle_tracer_mpas_ocean(MPI_COMM_WORLD);

#if NDARRAY_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
