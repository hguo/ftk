#include <fstream>
#include <mutex>
#include <set>
#include <cassert>
#include "ftk/external/cxxopts.hpp"
#include "ftk/ndarray/synthetic.hh"
#include "ftk/filters/parallel_vector_tracker_3d_regular.hh"
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
std::string archived_discrete_pvs_filenames, archived_traced_pvs_filename;
std::string accelerator;
int nthreads = std::thread::hardware_concurrency();
bool verbose = false, demo = false, show_vtk = false, help = false;

// input stream
ftk::ndarray_stream<> stream;

///////////////////////////////
int parse_arguments(int argc, char **argv)
{
  const int argc0 = argc;

  cxxopts::Options options(argv[0]);
  options.add_options()COMMON_OPTS_INPUTS()
    ("archived-discrete-parallel-vectors", "Archived discrete parallel vector points", cxxopts::value<std::string>(archived_discrete_pvs_filenames))
    ("archived-traced-parallel-vectors", "Archived traced critical points", cxxopts::value<std::string>(archived_traced_pvs_filename))
    ("o,output", "Output file, either one single file (e.g. out.vtp) or a pattern (e.g. out-%05d.vtp)", 
     cxxopts::value<std::string>(output_filename))
    ("output-type", "Output type {discrete|traced|sliced|intercepted}, by default traced", 
     cxxopts::value<std::string>(output_type)->default_value("traced"))
    ("output-format", "Output format {auto|text|vtp}, by default auto", 
     cxxopts::value<std::string>(output_format)->default_value(str_auto))
    ("nthreads", "Number of threads", 
     cxxopts::value<int>(nthreads))
    ("a,accelerator", "Accelerator {none|cuda} (experimental)",
     cxxopts::value<std::string>(accelerator)->default_value(str_none))
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

  if (output_filename.empty())
    fatal("Missing '--output'.");
  
  nlohmann::json j_input = args_to_json(results);
  stream.set_input_source_json(j_input);

  const auto js = stream.get_json();
  const size_t nd = stream.n_dimensions(),
               DW = js["dimensions"][0], 
               DH = js["dimensions"].size() > 1 ? js["dimensions"][1].get<int>() : 0,
               DD = js["dimensions"].size() > 2 ? js["dimensions"][2].get<int>() : 0;

  ftk::parallel_vector_tracker_3d_regular tracker;
  tracker.set_domain(ftk::lattice({1, 1, 1}, {DW-2, DH-2, DD-2}));
  tracker.set_array_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD}));
  tracker.initialize();

  stream.set_callback([&](int k, const ftk::ndarray<double> &field_data) {
    const double sigma = 0.6;


    tracker.push_field_data_snapshot(field_data, ftk::Jv_dot_v3(field_data));
    tracker.update_timestep();
  });
  stream.start();
  stream.finish();
  tracker.finalize();

  // tracker.write_discrete_pvs_vtk(output_filename);
  tracker.write_traced_pvs_vtk(output_filename);

  return 0;
}

int main(int argc, char **argv)
{
  diy::mpi::environment env;
  diy::mpi::communicator world;
  
  parse_arguments(argc, argv);

#if 0
  wrapper.consume(stream);
  if (!disable_post_processing)
     wrapper.post_process();
  wrapper.write();
#endif
 
  return 0;
}
