#include <fstream>
#include <mutex>
#include <set>
#include <cassert>
#include "ftk/external/cxxopts.hpp"
#include "ftk/ndarray/synthetic.hh"
#include "ftk/filters/critical_point_tracker_2d_regular.hh"
#include "ftk/filters/critical_point_tracker_3d_regular.hh"
#include "ftk/filters/streaming_filter.hh"
#include "ftk/ndarray.hh"
#include "ftk/ndarray/conv.hh"
#include "cli_constants.hh"

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
std::string output_filename,
  output_format;
std::string accelerator;
std::string type_filter_str;
int nthreads = std::thread::hardware_concurrency();
bool verbose = false, demo = false, show_vtk = false, help = false;
bool use_type_filter = false;
unsigned int type_filter = 0;

// tracker and input stream
ftk::critical_point_tracker_regular* tracker = NULL;
ftk::ndarray_stream<> stream;

// constants
static const std::set<std::string> set_valid_output_format({str_auto, str_text, str_vtp});

///////////////////////////////
int parse_arguments(int argc, char **argv)
{
  cxxopts::Options options(argv[0]);
  options.add_options()COMMON_OPTS_INPUTS()
    ("o,output", "Output file", 
     cxxopts::value<std::string>(output_filename))
    ("type-filter", "Type filter: ane single or a combination of critical point types, e.g. `min', `max', `saddle', `min|max'",
     cxxopts::value<std::string>(type_filter_str))
    ("r,output-format", "Output format (auto|text|vtp)", 
     cxxopts::value<std::string>(output_format)->default_value(str_auto))
    ("nthreads", "Number of threads", 
     cxxopts::value<int>(nthreads))
    ("a,accelerator", "Accelerator (none|cuda)",
     cxxopts::value<std::string>(accelerator)->default_value(str_none))
    ("vtk", "Show visualization with vtk", 
     cxxopts::value<bool>(show_vtk))
    ("v,verbose", "Verbose outputs", cxxopts::value<bool>(verbose))
    ("help", "Print usage", cxxopts::value<bool>(help));
  auto results = options.parse(argc, argv);

  auto fatal = [&](const std::string& str) {
    std::cerr << "FATAL: " << str << std::endl
              << options.help() << std::endl;
    exit(1);
  };
  auto warn = [&](const std::string& str) {
    std::cerr << "WARN: " << str << std::endl;
  };
  
  if (help) {
    std::cerr << options.help() << std::endl;
    exit(0); // return 0;
  }

  auto j = parse_input_json(results);
  stream.set_input_source_json(j);
 
  // sanity check of arguments
  if (set_valid_output_format.find(output_format) == set_valid_output_format.end())
    fatal("invalid '--output-format'");
  if (set_valid_accelerator.find(accelerator) == set_valid_accelerator.end())
    fatal("invalid '--accelerator'");

  if (type_filter_str.size() > 0) {
    if (type_filter_str.find(str_critical_point_type_min) != std::string::npos)
      type_filter |= ftk::CRITICAL_POINT_2D_MINIMUM;
    if (type_filter_str.find(str_critical_point_type_max) != std::string::npos)
      type_filter |= ftk::CRITICAL_POINT_2D_MAXIMUM;
    if (type_filter_str.find(str_critical_point_type_saddle) != std::string::npos)
      type_filter |= ftk::CRITICAL_POINT_2D_SADDLE;
    if (!type_filter) fatal("Invalid type filter");
    use_type_filter = true;
  }

  // processing output
  if (show_vtk) {
#if FTK_HAVE_VTK
#else
    fatal("FTK not compiled with VTK.");
#endif
  }
  
  if (!show_vtk) { // output is optional if results are visualized with vtk
    if (output_filename.empty())
      fatal("Missing '--output'.");
  }

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

  fprintf(stderr, "SUMMARY\n=============\n");
  std::cerr << "input=" << j << std::endl;
  fprintf(stderr, "output_filename=%s\n", output_filename.c_str());
  fprintf(stderr, "output_format=%s\n", output_format.c_str());
  fprintf(stderr, "type_filter=%s\n", type_filter_str.c_str());
  fprintf(stderr, "nthreads=%d\n", nthreads);
  fprintf(stderr, "=============\n");

  // assert(nd == 2 || nd == 3);
  // assert(nv == 1 || nv == 2 || nv == 3);
  // assert(DT > 0);

  return 0;
}

void track_critical_points()
{
  const auto j = stream.get_json();
  const size_t nd = j["nd"], DW = j["width"], DH = j["height"], DD = (nd == 2 ? 0 : j["depth"].get<size_t>()), DT = j["n_timesteps"];
  const size_t nv = stream.n_components();

  if (j["nd"] == 2) {
    tracker = new ftk::critical_point_tracker_2d_regular;
    tracker->set_array_domain(ftk::lattice({0, 0}, {DW, DH}));
  } else {
    tracker = new ftk::critical_point_tracker_3d_regular;
    tracker->set_array_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD}));
  }
  
  tracker->set_number_of_threads(nthreads);
      
  tracker->set_input_array_partial(false); // input data are not distributed

  if (use_type_filter)
    tracker->set_type_filter(type_filter);
  
  if (nv == 1) { // scalar field
    tracker->set_scalar_field_source( ftk::SOURCE_GIVEN );
    tracker->set_vector_field_source( ftk::SOURCE_DERIVED );
    tracker->set_jacobian_field_source( ftk::SOURCE_DERIVED );
    if (nd == 2) { // 2D
      tracker->set_domain(ftk::lattice({2, 2}, {DW-3, DH-3})); // the indentation is needed becase both gradient and jacoobian field will be automatically derived
    } else { // 3D
      tracker->set_domain(ftk::lattice({2, 2, 2}, {DW-3, DH-3, DD-3})); // the indentation is needed becase both gradient and jacoobian field will be automatically derived
    }
  } else { // vector field
    tracker->set_scalar_field_source( ftk::SOURCE_NONE );
    tracker->set_vector_field_source( ftk::SOURCE_GIVEN );
    tracker->set_jacobian_field_source( ftk::SOURCE_DERIVED );
    if (nd == 2) { // 2D
      tracker->set_domain(ftk::lattice({1, 1}, {DW-2, DH-2})); // the indentation is needed becase the jacoobian field will be automatically derived
    } else {
      tracker->set_domain(ftk::lattice({1, 1, 1}, {DW-2, DH-2, DD-2})); // the indentation is needed becase the jacoobian field will be automatically derived
    }
  }
  tracker->initialize();
   

  auto push_timestep = [&](const ftk::ndarray<double>& field_data) {
    if (nv == 1) { // scalar field
#if 0
      if (spatial_smoothing) {
        ftk::ndarray<double> scalar = 
          ftk::conv2D_gaussian(field_data, spatial_smoothing, 
              spatial_smoothing_kernel_size, spatial_smoothing_kernel_size, 2);
        tracker->push_scalar_field_snapshot(scalar);
      } else 
#endif
        tracker->push_scalar_field_snapshot(field_data);
    }
    else // vector field
      tracker->push_vector_field_snapshot(field_data);
  };

  stream.set_callback([&](int k, ftk::ndarray<double> field_data) {
    push_timestep(field_data);
    if (k != 0) tracker->advance_timestep();
    if (k == DT-1) {fprintf(stderr, "last!\n"); tracker->update_timestep();}
  });

  stream.start();
  stream.finish();
  tracker->finalize();
  // delete tracker;
}

void write_outputs()
{
  if (output_filename.empty()) return;

  if (output_format == str_vtp) 
    tracker->write_traced_critical_points_vtk(output_filename);
  else if (output_format == str_text) 
    tracker->write_traced_critical_points_text(output_filename);
}

void start_vtk_window()
{
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
  parse_arguments(argc, argv);
  track_critical_points();
    
  write_outputs();
  
  if (show_vtk) 
    start_vtk_window();

  delete tracker;
  return 0;
}
