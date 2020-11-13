#include <fstream>
#include <mutex>
#include <set>
#include <cassert>
#include "ftk/external/cxxopts.hpp"
#include "ftk/ndarray/synthetic.hh"
#include "ftk/filters/contour_tracker_2d_regular.hh"
#include "ftk/filters/contour_tracker_3d_regular.hh"
#include "ftk/filters/streaming_filter.hh"
#include "ftk/ndarray.hh"
#include "ftk/ndarray/conv.hh"
#include "constants.hh"

diy::mpi::environment env;
  
// global variables
std::string output_filename, output_type, output_format;
std::string archived_discrete_pvs_filenames, archived_traced_pvs_filename;
std::string accelerator;
int nthreads = std::thread::hardware_concurrency();
bool verbose = false, demo = false, show_vtk = false, help = false;

// determined later
double threshold = 0.0;

// input stream
ftk::ndarray_stream<> stream;

///////////////////////////////
int parse_arguments(int argc, char **argv)
{
  diy::mpi::communicator world;
  const int argc0 = argc;

  cxxopts::Options options(argv[0]);
  options.add_options()COMMON_OPTS_INPUTS()
    ("threshold", "Threshold for levelset tracking", cxxopts::value<double>(threshold))
    ("o,output", "Output file, either one single file (e.g. out.vtp) or a pattern (e.g. out-%05d.vtp)", 
     cxxopts::value<std::string>(output_filename))
    ("output-type", "Output type {isovolume|sliced|intersections}, by default isovolume", 
     cxxopts::value<std::string>(output_type)->default_value("isovolume"))
    ("output-format", "Output format {auto|text|vtp|vtu|ply}, by default auto", 
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
 
  if (world.rank() == 0) {
    fprintf(stderr, "SUMMARY\n=============\n");
    std::cerr << "input=" << stream.get_json() << std::endl;
    fprintf(stderr, "output_format=%s\n", output_format.c_str());
    fprintf(stderr, "threshold=%f\n", threshold);
    fprintf(stderr, "nthreads=%d\n", nthreads);
    fprintf(stderr, "=============\n");
  }

  const auto js = stream.get_json();
  const size_t nd = stream.n_dimensions(),
               DW = js["dimensions"][0], 
               DH = js["dimensions"].size() > 1 ? js["dimensions"][1].get<int>() : 0,
               DD = js["dimensions"].size() > 2 ? js["dimensions"][2].get<int>() : 0;
  const int nt = js["n_timesteps"];

  ftk::contour_tracker_regular *tracker;
  if (DD == 0) tracker = new ftk::contour_tracker_2d_regular(world);
  else tracker = new ftk::contour_tracker_3d_regular(world);
  
  if (accelerator == "cuda")
    tracker->use_accelerator(ftk::FTK_XL_CUDA);

  tracker->set_domain(ftk::lattice({0, 0, 0}, {DW-2, DH-2, DD-2}));
  tracker->set_array_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD}));
  tracker->set_end_timestep(nt - 1);
  tracker->set_number_of_threads(nthreads);
  tracker->set_threshold(threshold);
  tracker->initialize();

  stream.set_callback([&](int k, const ftk::ndarray<double> &field_data) {
    tracker->push_field_data_snapshot(field_data);
    
    if (k != 0) tracker->advance_timestep();
    if (k == nt - 1) tracker->update_timestep();
  });
  stream.start();
  stream.finish();
  tracker->finalize();

  if (output_type == "intersections") {
    tracker->write_intersections_vtp(output_filename);
  } else if (output_type == "sliced") {
    tracker->write_sliced_vtu(output_filename);
  } else {
    tracker->write_isovolume_vtu(output_filename);
  }
  
  delete tracker;
  return 0;
}

int main(int argc, char **argv)
{
  parse_arguments(argc, argv);
 
  return 0;
}
