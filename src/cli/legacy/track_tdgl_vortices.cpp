#include <fstream>
#include <mutex>
#include <set>
#include <cassert>
#include "ftk/external/cxxopts.hpp"
#include "ftk/ndarray/synthetic.hh"
#include "ftk/filters/tdgl_vortex_tracker_3d_regular.hh"
#include "ftk/filters/streaming_filter.hh"
#include "ftk/ndarray.hh"
#include "ftk/ndarray/conv.hh"
#include "constants.hh"
  
diy::mpi::environment env;

// global variables
std::string input_pattern;
std::string archived_intersections, archived_surfaces;
std::string output_pattern, output_type, output_format;
std::string accelerator;
size_t ntimesteps = 0;
int nthreads = std::thread::hardware_concurrency();
bool verbose = false, help = false;

// input stream
ftk::ndarray_stream<> stream;

///////////////////////////////
int parse_arguments(int argc, char **argv)
{
  const int argc0 = argc;

  cxxopts::Options options(argv[0]);
  options.add_options()
    ("i,input", "Input file name pattern", 
     cxxopts::value<std::string>(input_pattern))
    ("archived-intersections", "Archived vortex intersections", 
     cxxopts::value<std::string>(archived_intersections))
    ("archived-surfaces", "Archived vortex surfaces", 
     cxxopts::value<std::string>(archived_surfaces))
    ("n,timesteps", "Number of timesteps (override)", 
     cxxopts::value<size_t>(ntimesteps))
    ("o,output", "Output file, either one single file (e.g. out.vtp) or a pattern (e.g. out-%05d.vtp)", 
     cxxopts::value<std::string>(output_pattern))
    ("output-type", "Output type {surfaces|sliced|intersections}, by default surfaces", 
     cxxopts::value<std::string>(output_type)->default_value("surfaces"))
    ("output-format", "Output format {text|vtp|vtu|ply}",
     cxxopts::value<std::string>(output_format))
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

  if (output_pattern.empty())
    fatal("Missing '--output'.");
  
  auto filenames = ftk::ndarray<float>::glob(input_pattern);
  if (filenames.empty()) 
    fatal("unable to find matching filenames.");

  if (ntimesteps != 0) filenames.resize(ntimesteps);

  ftk::tdgl_reader meta_reader(filenames[0], false); // read metadata without read actuall data
  meta_reader.read();
  const auto &meta = meta_reader.get_meta();
  
  const size_t DW = meta.dims[0], DH = meta.dims[1], DD = meta.dims[2];
  const int nt = filenames.size();
 
  diy::mpi::communicator world;
  if (world.rank() == 0) {
    fprintf(stderr, "SUMMARY\n=============\n");
    // std::cerr << "input=" << stream.get_json() << std::endl;
    fprintf(stderr, "input_pattern=%s\n", input_pattern.c_str());
    fprintf(stderr, "dims=%zu, %zu, %zu\n", DW, DH, DD);
    fprintf(stderr, "nt=%d\n", nt);
    fprintf(stderr, "output_format=%s\n", output_format.c_str());
    fprintf(stderr, "nthreads=%d\n", nthreads);
    fprintf(stderr, "=============\n");
  }

  ftk::tdgl_vortex_tracker_3d_regular tracker(world);

  if (accelerator == "cuda")
    tracker.use_accelerator(ftk::FTK_XL_CUDA);

  tracker.set_domain(ftk::lattice({0, 0, 0}, {DW-2, DH-2, DD-2}));
  tracker.set_array_domain(ftk::lattice({0, 0, 0}, {DW, DH, DD}));
  tracker.set_end_timestep(nt - 1);
  tracker.set_number_of_threads(nthreads);
  tracker.initialize();

  if (!archived_intersections.empty()) {
    // TODO
  } else if (!archived_surfaces.empty()) {
    tracker.read_surfaces(archived_surfaces);
  } else {
    for (int k = 0; k < filenames.size(); k ++) {
      ftk::tdgl_reader reader(filenames[k]);
      reader.read();

      tracker.push_field_data_snapshot(reader.meta, 
          reader.rho, reader.phi, 
          reader.re, reader.im);
      
      if (k != 0) tracker.advance_timestep();
      if (k == nt - 1) tracker.update_timestep();
    }
    tracker.finalize();
  }

  if (output_type == "intersections")
    tracker.write_intersections(output_pattern);
  else if (output_type == "sliced")
    tracker.write_sliced(output_pattern);
  else 
    tracker.write_surfaces(output_pattern, output_format);
  
  return 0;
}

int main(int argc, char **argv)
{
  diy::mpi::communicator world;
  
  parse_arguments(argc, argv);
 
  return 0;
}
