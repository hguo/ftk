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

#if FTK_HAVE_VTK
#include <ftk/geometry/curve2vtk.hh>
#include <vtkPolyDataMapper.h>
#include <vtkTubeFilter.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#endif
  
// global variables
std::string output_filename, output_type, output_format;
std::string accelerator;
std::string type_filter_str;
int nthreads = std::thread::hardware_concurrency();
bool verbose = false, demo = false, show_vtk = false, help = false;
bool enable_streaming_trajectories = false, 
     enable_discarding_interval_points = false;

// tracker and input stream
ftk::critical_point_tracker_wrapper wrapper;
// ftk::critical_point_tracker_regular* tracker = NULL;
ftk::ndarray_stream<> stream;

nlohmann::json j_tracker;

// constants
static const std::set<std::string> set_valid_output_format({str_auto, str_text, str_vtp});

///////////////////////////////
int parse_arguments(int argc, char **argv)
{
  const int argc0 = argc;

  cxxopts::Options options(argv[0]);
  options.add_options()COMMON_OPTS_INPUTS()
    ("o,output", "Output file, either one single file (e.g. out.vtp) or a pattern (e.g. out-%05d.vtp)", 
     cxxopts::value<std::string>(output_filename))
    ("output-type", "Output type {discrete|traced|sliced}, by default traced", 
     cxxopts::value<std::string>(output_type)->default_value("traced"))
    ("output-format", "Output format {auto|text|vtp}, by default auto", 
     cxxopts::value<std::string>(output_format)->default_value(str_auto))
    ("type-filter", "Type filter: ane single or a combination of critical point types, e.g. `min', `max', `saddle', `min|max'",
     cxxopts::value<std::string>(type_filter_str))
    ("nthreads", "Number of threads", 
     cxxopts::value<int>(nthreads))
    ("a,accelerator", "Accelerator {none|cuda} (experimental)",
     cxxopts::value<std::string>(accelerator)->default_value(str_none))
    ("stream",  "Stream trajectories (experimental)",
     cxxopts::value<bool>(enable_streaming_trajectories))
    ("discard-interval-points", "Discard interval critical points (experimental)", 
     cxxopts::value<bool>(enable_discarding_interval_points))
    ("vtk", "Show visualization with vtk (legacy)", 
     cxxopts::value<bool>(show_vtk))
    ("v,verbose", "Verbose outputs", cxxopts::value<bool>(verbose))
    ("help", "Print usage", cxxopts::value<bool>(help));
  auto results = options.parse(argc, argv);

  if ((argc0 < 2) || help) {
    std::cerr << options.help() << std::endl;
    return 0;
  }

  auto j = parse_input_json(results);
  stream.set_input_source_json(j);

  auto fatal = [&](const std::string& str) {
	  std::cerr << "FATAL: " << str << std::endl
	            << options.help() << std::endl;
	  exit(1);
	};

	auto warn = [&](const std::string& str) {
	  std::cerr << "WARN: " << str << std::endl;
	};

  // sanity check of arguments
  if (set_valid_output_format.find(output_format) == set_valid_output_format.end())
    fatal("invalid '--output-format'");
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

  if (output_format == str_auto) {
    if (ends_with(output_filename, str_ext_vtp)) {
#if FTK_HAVE_VTK
      output_format = str_vtp;
#else
      fatal("FTK not compiled with VTK.");
#endif
    }
    else 
      output_format = str_text;
  }
  
  nlohmann::json jt;
  jt["output"] = output_filename;
  jt["output_type"] = output_type;
  jt["output_format"] = output_format;
  
  if (enable_streaming_trajectories) 
    jt["enable_streaming_trajectories"] = true;

  if (enable_discarding_interval_points)
    jt["enable_discarding_interval_points"] = true;

  jt["type_filter"] = type_filter_str;

  wrapper.configure(jt);

  diy::mpi::communicator world;
  if (world.rank() == 0) {
    fprintf(stderr, "SUMMARY\n=============\n");
    std::cerr << "input=" << stream.get_json() << std::endl;
    std::cerr << "config=" << wrapper.get_json() << std::endl;
    // fprintf(stderr, "output_filename=%s\n", output_filename.c_str());
    // fprintf(stderr, "output_format=%s\n", output_format.c_str());
    // fprintf(stderr, "type_filter=%s\n", type_filter_str.c_str());
    fprintf(stderr, "nthreads=%d\n", nthreads);
    fprintf(stderr, "=============\n");
  }

  // assert(nd == 2 || nd == 3);
  // assert(nv == 1 || nv == 2 || nv == 3);
  // assert(DT > 0);

  return 0;
}

void start_vtk_window()
{
  auto tracker = wrapper.get_tracker();
#if FTK_HAVE_VTK
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
  diy::mpi::environment env;
  diy::mpi::communicator world;
  
  parse_arguments(argc, argv);
  wrapper.consume(stream);
 
  if (world.rank() == 0 && show_vtk) 
    start_vtk_window();

  return 0;
}