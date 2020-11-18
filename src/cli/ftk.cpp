#include <fstream>
#include <mutex>
#include <set>
#include <cassert>
#include "ftk/external/cxxopts.hpp"
#include "ftk/ndarray/synthetic.hh"
#include "ftk/filters/critical_point_tracker_2d_regular.hh"
#include "ftk/filters/critical_point_tracker_3d_regular.hh"
#include "ftk/filters/critical_point_tracker_wrapper.hh"
#include "ftk/filters/contour_tracker_2d_regular.hh"
#include "ftk/filters/contour_tracker_3d_regular.hh"
#include "ftk/filters/parallel_vector_tracker_3d_regular.hh"
#include "ftk/filters/tdgl_vortex_tracker_3d_regular.hh"
#include "ftk/filters/threshold_tracker.hh"
#include "ftk/filters/streaming_filter.hh"
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
std::string accelerator;
std::string type_filter_str;
int nthreads = std::thread::hardware_concurrency();
bool verbose = false, timing = false, help = false;
int nblocks; 
bool enable_streaming_trajectories = false, 
     enable_discarding_interval_points = false,
     enable_deriving_velocities = false,
     disable_robust_detection = false,
     disable_post_processing = false;
int intercept_length = 2;
double duration_pruning_threshold = 0.0;

size_t ntimesteps = 0;
std::vector<std::string> input_filenames;

// contour/levelset specific
double threshold = 0.0;

// xgc specific
std::string xgc_mesh_filename, 
  xgc_smoothing_kernel_filename = "xgc.kernel",
  xgc_write_back_filename;
bool xgc_post_process = false, 
     xgc_torus = false;
double xgc_smoothing_kernel_size = 0.03;

// tracker and input stream
std::shared_ptr<critical_point_tracker_wrapper> tracker_critical_point;
std::shared_ptr<contour_tracker_regular> tracker_contour;
std::shared_ptr<tdgl_vortex_tracker_3d_regular> tracker_tdgl;
std::shared_ptr<threshold_tracker<>> tracker_threshold;
std::shared_ptr<ndarray_stream<>> stream;

nlohmann::json j_input, j_tracker;

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
        str_scalar("scalar"),
        str_vector("vector"),
        str_text("text"),
        str_cuda("cuda");

static const std::string
        str_ext_vti(".vti"), // vtkImageData
        str_ext_vtp(".vtp"), // vtkPolyData
        str_ext_ply(".ply"),
        str_ext_stl(".stl"),
        str_ext_netcdf(".nc"),
        str_ext_hdf5(".h5"),
        str_ext_adios2(".bp");

static const std::string
        str_critical_point_type_min("min"),
        str_critical_point_type_max("max"),
        str_critical_point_type_saddle("saddle");

static const std::string
        str_feature_critical_point("cp"),
        str_feature_tdgl_vortex("tdgl"),
        str_feature_isosurface("iso");

static const std::set<std::string>
        set_valid_accelerator({str_none, str_cuda}),
        set_valid_input_format({str_auto, str_float32, str_float64, str_netcdf, str_hdf5, str_vti, str_adios2}),
        set_valid_input_dimension({str_auto, str_two, str_three});

static void fatal(const cxxopts::Options &options, const std::string& str) {
  std::cerr << "FATAL: " << str << std::endl
            << options.help() << std::endl;
  exit(1);
};

static void fatal(const std::string& str) {
  std::cerr << "FATAL: " << str << std::endl;
  exit(1);
};

void warn(const std::string& str) {
  std::cerr << "WARN: " << str << std::endl;
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

  if (accelerator != str_none)
    j_tracker["accelerator"] = accelerator;

  if (archived_intersections_filename.size() > 0)
    j_tracker["archived_discrete_critical_points_filename"] = archived_intersections_filename;
  
  if (archived_traced_filename.size() > 0)
    j_tracker["archived_traced_critical_points_filename"] = archived_traced_filename;
  
  if (enable_streaming_trajectories) 
    j_tracker["enable_streaming_trajectories"] = true;

  if (enable_discarding_interval_points)
    j_tracker["enable_discarding_interval_points"] = true;

  if (enable_deriving_velocities)
    j_tracker["enable_deriving_velocities"] = true;

  if (disable_robust_detection)
    j_tracker["enable_robust_detection"] = false;

  if (duration_pruning_threshold > 0)
    j_tracker["duration_pruning_threshold"] = duration_pruning_threshold;

  j_tracker["enable_post_processing"] = !disable_post_processing;

  j_tracker["type_filter"] = type_filter_str;

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

  tracker_critical_point.reset(new critical_point_tracker_wrapper);
  tracker_critical_point->configure(j_tracker);

  if (comm.rank() == 0) {
    // fprintf(stderr, "SUMMARY\n=============\n");
    std::cerr << "input=" << std::setw(2) << stream->get_json() << std::endl;
    std::cerr << "config=" << std::setw(2) << tracker_critical_point->get_json() << std::endl;
    // fprintf(stderr, "=============\n");
  }

  // assert(nd == 2 || nd == 3);
  // assert(nv == 1 || nv == 2 || nv == 3);
  // assert(DT > 0);
}

static void execute_critical_point_tracker(diy::mpi::communicator comm)
{
  tracker_critical_point->consume(*stream);
 
  if (!disable_post_processing)
     tracker_critical_point->post_process();
  
  tracker_critical_point->write();
}

void initialize_contour_tracker(diy::mpi::communicator comm)
{
  const auto js = stream->get_json();
  const size_t nd = stream->n_dimensions(),
               DW = js["dimensions"][0], 
               DH = js["dimensions"].size() > 1 ? js["dimensions"][1].get<int>() : 0,
               DD = js["dimensions"].size() > 2 ? js["dimensions"][2].get<int>() : 0;
  const int nt = js["n_timesteps"];

  if (DD == 0) tracker_contour.reset( new contour_tracker_2d_regular(comm) ); 
  else tracker_contour.reset(new contour_tracker_3d_regular(comm) );
  
  if (accelerator == "cuda")
    tracker_contour->use_accelerator(FTK_XL_CUDA);

  tracker_contour->set_domain(lattice({0, 0, 0}, {DW-2, DH-2, DD-2}));
  tracker_contour->set_array_domain(lattice({0, 0, 0}, {DW, DH, DD}));
  tracker_contour->set_end_timestep(nt - 1);
  tracker_contour->set_number_of_threads(nthreads);
  tracker_contour->set_threshold(threshold);
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

void initialize_tdgl_tracker(diy::mpi::communicator comm)
{
  input_filenames = ndarray<float>::glob(input_pattern);
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

  if (results.count("synthetic")) {
    j["type"] = "synthetic";
    j["name"] = results["synthetic"].as<std::string>();
  } else 
    j["type"] = "file";

  if (results.count("input-format")) j["format"] = results["input-format"].as<std::string>();
  if (results.count("dim")) j["nd"] = results["dim"].as<std::string>();
  
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
    ("n,timesteps", "Number of timesteps", cxxopts::value<size_t>(ntimesteps))
    ("var", "Variable name(s), e.g. `scalar', `u,v,w'.  Valid only for NetCDF, HDF5, and VTK.", cxxopts::value<std::string>())
    ("temporal-smoothing-kernel", "Temporal smoothing kernel bandwidth", cxxopts::value<double>())
    ("temporal-smoothing-kernel-size", "Temporal smoothing kernel size", cxxopts::value<size_t>())
    ("spatial-smoothing-kernel", "Spatial smoothing kernel bandwidth", cxxopts::value<double>())
    ("spatial-smoothing-kernel-size", "Spatial smoothing kernel size", cxxopts::value<size_t>())
    ("perturbation", "Gaussian perturbation sigma", cxxopts::value<double>())
    ("m,mesh", "Input mesh file (will shadow arguments including width, height, depth)", cxxopts::value<std::string>())
    ("nblocks", "Number of total blocks", cxxopts::value<int>(nblocks))
    // ("archived-discrete-critical-points", "Archived discrete critical points", cxxopts::value<std::string>(archived_discrete_critical_points_filename))
    ("archived-intersections", "Archived discrete intersections", cxxopts::value<std::string>(archived_intersections_filename))
    ("archived-traced", "Archived traced results", cxxopts::value<std::string>(archived_traced_filename))
    ("xgc-mesh", "XGC mesh file", cxxopts::value<std::string>(xgc_mesh_filename))
    ("xgc-smoothing-kernel-file", "XGC: smoothing kernel file", cxxopts::value<std::string>(xgc_smoothing_kernel_filename))
    ("xgc-smoothing-kernel-size", "XGC: smoothing kernel size", cxxopts::value<double>(xgc_smoothing_kernel_size))
    ("xgc-torus", "XGC: track over poloidal planes", cxxopts::value<bool>(xgc_torus))
    ("xgc-write-back", "XGC: write original back into vtu files", cxxopts::value<std::string>(xgc_write_back_filename))
    ("xgc-post-process", "XGC: enable post-processing", cxxopts::value<bool>(xgc_post_process))
    ("o,output", "Output file, either one single file (e.g. out.vtp) or a pattern (e.g. out-%05d.vtp)", 
     cxxopts::value<std::string>(output_pattern))
    ("output-type", "Output type {discrete|traced|sliced|intercepted}, by default traced", 
     cxxopts::value<std::string>(output_type)->default_value("traced"))
    ("output-format", "Output format {text|vtp|vtu|ply}.  The default behavior is to automatically determine format by filename", 
     cxxopts::value<std::string>(output_format)->default_value(str_auto))
    ("intercept-length", "Length of intercepted outputs", 
     cxxopts::value<int>(intercept_length)->default_value("2"))
    ("type-filter", "Type filter: ane single or a combination of critical point types, e.g. `min', `max', `saddle', `min|max'",
     cxxopts::value<std::string>(type_filter_str))
    ("nthreads", "Number of threads", 
     cxxopts::value<int>(nthreads))
    ("timing", "Enable timing", 
     cxxopts::value<bool>(timing))
    ("a,accelerator", "Accelerator {none|cuda} (experimental)",
     cxxopts::value<std::string>(accelerator)->default_value(str_none))
    ("stream",  "Stream trajectories (experimental)",
     cxxopts::value<bool>(enable_streaming_trajectories))
    ("discard-interval-points", "Discard interval critical points (experimental)", 
     cxxopts::value<bool>(enable_discarding_interval_points))
    ("derive-velocities", "Enable deriving velocities", 
     cxxopts::value<bool>(enable_deriving_velocities))
    ("no-robust-detection", "Disable robust detection (faster than robust detection)",
     cxxopts::value<bool>(disable_robust_detection))
    ("no-post-processing", "Disable post-processing",
     cxxopts::value<bool>(disable_post_processing))
    ("duration-pruning", "Prune trajectories below certain duration", 
     cxxopts::value<double>(duration_pruning_threshold))
    ("v,verbose", "Verbose outputs", cxxopts::value<bool>(verbose))
    ("help", "Print usage", cxxopts::value<bool>(help));
  auto results = options.parse(argc, argv);

  // sanity check of arguments
  if (set_valid_accelerator.find(accelerator) == set_valid_accelerator.end())
    fatal(options, "invalid '--accelerator'");

  if (output_pattern.empty())
    fatal(options, "Missing '--output'.");
  
  ttype = tracker::str2tracker(feature);
  if (ttype != TRACKER_TDGL_VORTEX) { // TDGL uses a different reader for now
    j_input = args_to_input_stream_json(results);
    stream->set_input_source_json(j_input);
  }

  if (ttype == TRACKER_CRITICAL_POINT)
    initialize_critical_point_tracker(comm);
  else if (ttype == TRACKER_CONTOUR)
    initialize_contour_tracker(comm); 
  else if (ttype == TRACKER_TDGL_VORTEX)
    initialize_tdgl_tracker(comm);
  else 
    fatal(options, "missing or invalid '--feature'");

  return 0;
}

int main(int argc, char **argv)
{
  diy::mpi::environment env(argc, argv);
  diy::mpi::communicator comm;
  
  stream.reset(new ndarray_stream<>);
 
  parse_arguments(argc, argv, comm);

  if (ttype == TRACKER_CRITICAL_POINT)
    execute_critical_point_tracker(comm);
  else if (ttype == TRACKER_CONTOUR)
    execute_contour_tracker(comm);
  else if (ttype == TRACKER_TDGL_VORTEX)
    execute_tdgl_tracker(comm);

  return 0;
}
