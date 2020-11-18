#include <fstream>
#include <mutex>
#include <set>
#include <cassert>
#include <iostream>
#include <ftk/external/cxxopts.hpp>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/filters/threshold_tracker.hh>
#include <ftk/filters/connected_component_tracker.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/stream.hh>
#include <ftk/algorithms/hoshen_kopelman.hh>
#include <ftk/tracking_graph/tracking_graph.hh>
#include "constants.hh"

diy::mpi::environment env;
diy::mpi::communicator world;

using nlohmann::json;

std::string output_filename_pattern,
  output_format,
  output_filename_dot;
std::string accelerator;
int nthreads = std::thread::hardware_concurrency();
bool verbose = false, help = false;

// determined later
double threshold = 0.0;

// constants
static const std::set<std::string> set_valid_output_format({str_auto, str_text, str_vti});

// tracker
ftk::threshold_tracker<> tracker(world);
ftk::ndarray_stream<> stream;

int parse_arguments(int argc, char **argv)
{
  cxxopts::Options options(argv[0]);
  options.add_options()COMMON_OPTS_INPUTS()
    ("threshold", "Threshold for threshold tracking", cxxopts::value<double>(threshold))
    ("o,output", "Output file name pattern, e.g. 'out-%d.raw', 'out-%04d.vti'",
     cxxopts::value<std::string>(output_filename_pattern))
    ("output-format", "Output file format (auto|raw|nc|h5|vti)",
     cxxopts::value<std::string>(output_format)->default_value(str_auto))
    ("write-graph-dot", "Write tracking graph in GraphViz format",
     cxxopts::value<std::string>(output_filename_dot))
    ("nthreads", "Number of threads", 
     cxxopts::value<int>(nthreads))
    ("a,accelerator", "Accelerator (none|cuda)",
     cxxopts::value<std::string>(accelerator)->default_value(str_none))
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

  // sanity check of arguments
  if (set_valid_output_format.find(output_format) == set_valid_output_format.end())
    fatal("invalid '--output-format'");
  if (set_valid_accelerator.find(accelerator) == set_valid_accelerator.end())
    fatal("invalid '--accelerator'");
 
  if (!results.count("threshold")) 
    fatal("Missing threshold value.");

  // configure input stream
  json j_input = args_to_json(results);
  stream.set_input_source_json(j_input);

  if (output_filename_pattern.empty())
    fatal("Missing '--output'.");

  if (output_format == str_auto) {
    if (ftk::ends_with(output_filename_pattern, str_ext_vtp)) {
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
  std::cerr << "input=" << stream.get_json() << std::endl;
  fprintf(stderr, "output_filename_pattern=%s\n", output_filename_pattern.c_str());
  fprintf(stderr, "output_format=%s\n", output_format.c_str());
  fprintf(stderr, "threshold=%f\n", threshold);
  fprintf(stderr, "nthreads=%d\n", nthreads);
  fprintf(stderr, "=============\n");

  // assert(nd == 2 || nd == 3);
  // assert(DT > 0);

  return 0;
}

void track_threshold()
{
  tracker.set_threshold( threshold );

  stream.set_callback([&](int k, ftk::ndarray<double> field_data) {
    fprintf(stderr, "current_timestep=%d\n", k);
    tracker.push_scalar_field_data_snapshot(field_data);
    tracker.advance_timestep();
  });

  stream.start();
  stream.finish();
  tracker.finalize();

  const auto &tg = tracker.get_tracking_graph();
  tg.generate_dot_file("dot");
}

void write_outputs()
{
  // TODO
}

int main(int argc, char **argv)
{
  diy::mpi::environment env;
  parse_arguments(argc, argv);
  track_threshold();
  return 0;
}
