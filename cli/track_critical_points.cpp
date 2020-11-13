#include <fstream>
#include <mutex>
#include <set>
#include <cassert>
#include "ftk/external/cxxopts.hpp"
#include "ftk/ndarray/synthetic.hh"
#include "ftk/filters/critical_point_tracker_2d_regular.hh"
#include "ftk/filters/critical_point_tracker_3d_regular.hh"
#include "ftk/filters/critical_point_tracker_wrapper.hh"
#include "ftk/filters/streaming_filter.hh"
#include "ftk/ndarray.hh"
#include "ftk/ndarray/conv.hh"
#include "constants.hh"

#if FTK_HAVE_VTK && !FTK_HAVE_PARAVIEW
#include <ftk/geometry/curve2vtk.hh>
#include <vtkPolyDataMapper.h>
#include <vtkTubeFilter.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#endif
  
diy::mpi::environment env;
  
// global variables
std::string output_filename, output_type, output_format;
std::string mesh_filename;
std::string archived_discrete_critical_points_filename,
  archived_traced_critical_points_filename;
std::string accelerator;
std::string type_filter_str;
int nthreads = std::thread::hardware_concurrency();
bool verbose = false, timing = false, show_vtk = false, help = false;
int nblocks; 
bool enable_streaming_trajectories = false, 
     enable_discarding_interval_points = false,
     enable_deriving_velocities = false,
     disable_robust_detection = false,
     disable_post_processing = false;
int intercept_length = 2;
double duration_pruning_threshold = 0.0;

// xgc specific
std::string xgc_mesh_filename, 
  xgc_smoothing_kernel_filename = "xgc.kernel",
  xgc_write_back_filename;
bool xgc_post_process = false, 
     xgc_torus = false;
double xgc_smoothing_kernel_size = 0.03;

// tracker and input stream
ftk::critical_point_tracker_wrapper wrapper;
// ftk::critical_point_tracker_regular* tracker = NULL;
ftk::ndarray_stream<> stream;

nlohmann::json j_tracker;

///////////////////////////////
int parse_arguments(int argc, char **argv)
{
  diy::mpi::communicator world;
  const int argc0 = argc;

  cxxopts::Options options(argv[0]);
  options.add_options()COMMON_OPTS_INPUTS()
    ("m,mesh", "Input mesh file (will shadow arguments including width, height, depth)", cxxopts::value<std::string>())
    ("nblocks", "Number of total blocks", cxxopts::value<int>(nblocks))
    ("archived-discrete-critical-points", "Archived discrete critical points", cxxopts::value<std::string>(archived_discrete_critical_points_filename))
    ("archived-traced-critical-points", "Archived discrete critical points", cxxopts::value<std::string>(archived_traced_critical_points_filename))
    ("xgc-mesh", "XGC mesh file", cxxopts::value<std::string>(xgc_mesh_filename))
    ("xgc-smoothing-kernel-file", "XGC: smoothing kernel file", cxxopts::value<std::string>(xgc_smoothing_kernel_filename))
    ("xgc-smoothing-kernel-size", "XGC: smoothing kernel size", cxxopts::value<double>(xgc_smoothing_kernel_size))
    ("xgc-torus", "XGC: track over poloidal planes", cxxopts::value<bool>(xgc_torus))
    ("xgc-write-back", "XGC: write original back into vtu files", cxxopts::value<std::string>(xgc_write_back_filename))
    ("xgc-post-process", "XGC: enable post-processing", cxxopts::value<bool>(xgc_post_process))
    ("o,output", "Output file, either one single file (e.g. out.vtp) or a pattern (e.g. out-%05d.vtp)", 
     cxxopts::value<std::string>(output_filename))
    ("output-type", "Output type {discrete|traced|sliced|intercepted}, by default traced", 
     cxxopts::value<std::string>(output_type)->default_value("traced"))
    ("output-format", "Output format {auto|text|vtp}, by default auto", 
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
    ("vtk", "Show visualization with vtk (legacy)", 
     cxxopts::value<bool>(show_vtk))
    ("v,verbose", "Verbose outputs", cxxopts::value<bool>(verbose))
    ("help", "Print usage", cxxopts::value<bool>(help));
  auto results = options.parse(argc, argv);

  if ((argc0 < 2) || help) {
    std::cerr << options.help() << std::endl;
    return 0;
  }

  auto fatal = [&](const std::string& str) {
	  std::cerr << "FATAL: " << str << std::endl
	            << options.help() << std::endl;
	  exit(1);
	};

	auto warn = [&](const std::string& str) {
	  std::cerr << "WARN: " << str << std::endl;
	};

  // sanity check of arguments
  if (set_valid_accelerator.find(accelerator) == set_valid_accelerator.end())
    fatal("invalid '--accelerator'");

  // processing output
  if (show_vtk) {
#if FTK_HAVE_VTK
#else
    fatal("FTK not compiled with VTK.");
#endif
  }

  if (output_filename.empty())
    fatal("Missing '--output'.");

  nlohmann::json j_input = args_to_json(results), 
                 j_tracker;
  
  stream.set_input_source_json(j_input);
  
  j_tracker["output"] = output_filename;
  j_tracker["output_type"] = output_type;
  j_tracker["output_format"] = output_format;
  j_tracker["intercept_length"] = intercept_length;

  j_tracker["nthreads"] = nthreads;
  j_tracker["enable_timing"] = timing;

  j_tracker["nblocks"] = std::max(world.size(), nblocks);

  if (accelerator != str_none)
    j_tracker["accelerator"] = accelerator;

  if (archived_discrete_critical_points_filename.size() > 0)
    j_tracker["archived_discrete_critical_points_filename"] = archived_discrete_critical_points_filename;
  
  if (archived_traced_critical_points_filename.size() > 0)
    j_tracker["archived_traced_critical_points_filename"] = archived_traced_critical_points_filename;
  
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

  wrapper.configure(j_tracker);

  if (world.rank() == 0) {
    // fprintf(stderr, "SUMMARY\n=============\n");
    std::cerr << "input=" << std::setw(2) << stream.get_json() << std::endl;
    std::cerr << "config=" << std::setw(2) << wrapper.get_json() << std::endl;
    // fprintf(stderr, "=============\n");
  }

  // assert(nd == 2 || nd == 3);
  // assert(nv == 1 || nv == 2 || nv == 3);
  // assert(DT > 0);

  return 0;
}

void start_vtk_window()
{
#if FTK_HAVE_VTK && !FTK_HAVE_PARAVIEW
  auto tracker = wrapper.get_tracker();
  auto vtkcurves = tracker->get_traced_critical_points_vtk();
  vtkcurves->Print(std::cerr);

  vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
  tubeFilter->SetInputData(vtkcurves);
  tubeFilter->SetRadius(1);
  tubeFilter->SetNumberOfSides(50);
  tubeFilter->Update();
  
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  // mapper->SetInputData(vtkcurves);
  mapper->SetInputConnection(tubeFilter->GetOutputPort());

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  // a renderer and render window
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  // add the actors to the scene
  renderer->AddActor(actor);
  renderer->SetBackground(1, 1, 1); // Background color white

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = 
      vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle( style );

  renderWindowInteractor->Start();
#else
  assert(false);
#endif
}

int main(int argc, char **argv)
{
  diy::mpi::communicator world;
  
  parse_arguments(argc, argv);
  wrapper.consume(stream);
 
  if (!disable_post_processing)
     wrapper.post_process();
  
  wrapper.write();
 
  if (world.rank() == 0 && show_vtk) 
    start_vtk_window();

  return 0;
}
