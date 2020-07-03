#include <fstream>
#include <mutex>
#include <set>
#include <cassert>
#include <iostream>
#include <ftk/external/cxxopts.hpp>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/filters/levelset_tracker.hh>
#include <ftk/filters/connected_component_tracker.hh>
#include <ftk/ndarray.hh>
#include <ftk/io/data_stream.hh>
#include <ftk/algorithms/hoshen_kopelman.hh>
#include <ftk/tracking_graph/tracking_graph.hh>
#include "cli_constants.hh"

std::string input_filename_pattern, 
  input_format,
  input_dimension,
  input_variable_name;
std::string output_filename_pattern,
  output_format,
  output_filename_dot;
std::vector<std::string> input_filenames; // assuming each file contains only one timestep, and all files have the exactly same structure
std::string accelerator;
size_t DW = 0, DH = 0, DD = 0, DT = 0;
int nthreads = std::thread::hardware_concurrency();
bool verbose = false, demo = false, show_vtk = false, help = false;

// determined later
int nd; // dimensionality
int varid = -1; // only for vti and netcdf
int ncdims = 0;  // Only for netcdf.  
                 // ndndims equals to either nd or (nd+1). 
                 // In the latter case, one of the netcdf dimension 
                 // is time.
int dimids[4] = {-1}; // Only for netcdf
size_t dimlens[4] = {0}; // Only for netcdf
double threshold = 0.0;

// constants
static const std::set<std::string> set_valid_output_format({str_auto, str_text, str_vti});

// tracker
ftk::levelset_tracker<> tracker;

ftk::ndarray<double> request_timestep(int k) // requesting k-th timestep
{
  // const double t = DT == 1 ? 0.0 : double(k)/(DT-1);
  // return ftk::synthetic_woven_2D<double>(DW, DH, t);
  const double t = DT == 1 ? 0.0 : double(k)/(DT-1) * 10;
  return ftk::synthetic_merger_2D<double>(DW, DH, t);
}

int parse_arguments(int argc, char **argv)
{
  cxxopts::Options options(argv[0]);
  options.add_options()
    ("i,input", "Input file name pattern: a single file or a series of file, e.g. 'scalar.raw', 'cm1out_000*.nc'", 
     cxxopts::value<std::string>(input_filename_pattern))
    ("f,input-format", "Input file format (auto|float32|float64|nc|h5|vti)", 
     cxxopts::value<std::string>(input_format)->default_value(str_auto))
    ("demo", "Use synthetic data in lieu of input files", cxxopts::value<bool>(demo))
    ("dim", "Spatial dimensionality of data (auto|2|3)", 
     cxxopts::value<std::string>(input_dimension)->default_value(str_auto))
    ("w,width", "Width", cxxopts::value<size_t>(DW))
    ("h,height", "Height", cxxopts::value<size_t>(DH))
    ("d,depth", "Depth", cxxopts::value<size_t>(DH))
    ("n,timesteps", "Number of timesteps", cxxopts::value<size_t>(DT))
    ("var", "Variable name (only for NetCDF, HDF5, and VTK)", 
     cxxopts::value<std::string>(input_variable_name))
    ("threshold", "Threshold for levelset tracking", 
     cxxopts::value<double>(threshold))
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
 
  // sanity check of arguments
  if (set_valid_input_format.find(input_format) == set_valid_input_format.end())
    fatal("invalid '--input-format'");
  if (set_valid_output_format.find(output_format) == set_valid_output_format.end())
    fatal("invalid '--output-format'");
  if (set_valid_accelerator.find(accelerator) == set_valid_accelerator.end())
    fatal("invalid '--accelerator'");
 
  if (input_dimension == str_auto || input_dimension.size() == 0) nd = 0; // auto
  else if (input_dimension == str_two) nd = 2;
  else if (input_dimension == str_three) nd = 3;
  else fatal("Unsupported data dimensionality.");

  if (!results.count("threshold")) 
    fatal("Missing threshold value.");

  // processing output
  if (show_vtk) {
#if FTK_HAVE_VTK
#else
    fatal("FTK not compiled with VTK.");
#endif
  }
  
  if (!show_vtk) { // output is optional if results are visualized with vtk
    if (output_filename_pattern.empty())
      fatal("Missing '--output'.");
  }

  if (output_format == str_auto) {
    if (ends_with(output_filename_pattern, str_ext_vtp)) {
#if FTK_HAVE_VTK
      output_format = str_vtp;
#else
      fatal("FTK not compiled with VTK.");
#endif
    }
    else 
      output_format = str_text;
  }

  // processing input
  if (demo) {
    if (input_filename_pattern.size())
      fprintf(stderr, "WARN: '--input' ignored.\n");
    if (input_format != "auto")
      fprintf(stderr, "WARN: '--input-format' ignored.\n");

    if (input_dimension == "3") nd = 3;
    else nd = 2;

    if (input_variable_name.size())
      warn("--var ignored.");

    // default DW, DH, DD, DT
    if (DW == 0) DW = 32;
    if (DH == 0) DH = 32;
    if (nd == 3 && DD == 0) DD = 32;
    if (DT == 0) DT = 32;
  } else {
    if (input_filename_pattern.size() == 0)
      fatal("'--input' empty.  Please specify inputs or use '--demo'");

    input_filenames = ftk::ndarray<double>::glob(input_filename_pattern);
    if (input_filenames.size() == 0)
      fatal("FATAL: cannot find input files");

    // determine input format
    if (input_format == str_auto) {
      if (ends_with(input_filenames[0], str_ext_vti))
        input_format = str_vti;
      else if (ends_with(input_filenames[0], str_ext_netcdf)) 
        input_format = str_netcdf;
      else if (ends_with(input_filenames[0], str_ext_hdf5)) 
        input_format = str_hdf5;
      else 
        fatal("Unable to determine file format.");
    }

    if (input_format == str_float32 || input_format == str_float64) {
      // determine spatial dimensionality 
      if (input_dimension == str_two) nd = 2;
      else if (input_dimension == str_three) nd = 3;
      else 
        fatal("Unable to determine spatial dimensionality.");

      // determine dimensions
      if ((nd == 2 && (DW == 0 || DH == 0)) ||
          (nd == 3 && (DW == 0 || DH == 0 || DW == 0)))
          fatal("Data dimensions are not specified.");

      // ignore var names
      if (input_variable_name.size())
        warn("Ignoring variable names.");
    } else if (input_format == str_vti) {
#if FTK_HAVE_VTK
      vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
      reader->SetFileName(input_filenames[0].c_str());
      reader->Update();

      vtkSmartPointer<vtkImageData> image = reader->GetOutput();

      // determine dimensionality
      int imageNd = image->GetDataDimension();
      if (input_dimension == str_auto) {
        nd = imageNd;
        if (!(nd == 2 || nd == 3))
          fatal("Unsupported VTI data dimensionality.");
      } else {
        if (nd != imageNd)
          fatal("Data dimensionality and given dimensionality mismatch.");
      }

      // determine dimensions
      if (DW || DH || DD)
        warn("Given data dimensions are ignored.");

      DW = image->GetDimensions()[0];
      DH = image->GetDimensions()[1];
      if (nd == 3)
        DD = image->GetDimensions()[2];

      image->GetPointData()->GetArray(input_variable_name.c_str(), varid);

      const int nv = image->GetPointData()->GetArray(varid)->GetNumberOfComponents();
      if (nv != 1)
        fatal("Number of components is not 1 in VTI.");

      // determine DT
      if (DT == 0) DT = input_filenames.size();
      else DT = std::min(DT, input_filenames.size());
#else
      fatal("FTK not compiled with VTK.");
#endif
    } else if (input_format == str_netcdf) {
#if FTK_HAVE_NETCDF
      if (input_variable_name.size() == 0)
        fatal("Variable name missing for NetCDF files.");

      int ncid;
      NC_SAFE_CALL( nc_open(input_filenames[0].c_str(), NC_NOWRITE, &ncid) );
      NC_SAFE_CALL( nc_inq_varid(ncid, input_variable_name.c_str(), &varid) );
      NC_SAFE_CALL( nc_inq_varndims(ncid, varid, &ncdims) );
      
      // determine dimensionality (nd)
      if (nd == 0) { // auto
        if (ncdims == 2) nd = 2;
        else if (ncdims == 3)
          fatal("Please specify the number of spatial dimensions with '--dim'");
        else if (ncdims == 4) 
          nd = 3;
        else 
          fatal("Unsupported NetCDF data dimensionality.");
      } else if (nd == ncdims) { // netcdf dimensions are spatial only
      } else if (nd == ncdims - 1) { // netcdf file has time dimension
        // NOTE: we assume the time dimension is NC_UNLIMITED, and each
        //       NetCDF contains only one time step.  Will remove this 
        //       limitation in future versions.
      } else {
        fprintf(stderr, "nd=%d, ncdims=%d\n", nd, ncdims);
        fatal("Unsupported NetCDF variable dimensionality.");
      }

      // determine spatial dimensions
      NC_SAFE_CALL( nc_inq_vardimid(ncid, varid, dimids) );
      for (int i = 0; i < ncdims; i ++) {
        NC_SAFE_CALL( nc_inq_dimlen(ncid, dimids[i], &dimlens[i]) );
        // fprintf(stderr, "dimlens %d: %d\n", i, dimlens[i]);
      }

      if (ncdims == 4) { // 3 spatial dims + 1 time dimension
        DW = dimlens[3];
        DH = dimlens[2];
        DD = dimlens[1];
      } else if (ncdims == 3) {
        if (nd == 3) {
          DW = dimlens[2];
          DH = dimlens[1];
          DD = dimlens[0];
        } else { // nd == 2
          DW = dimlens[2];
          DH = dimlens[1];
          DD = 0;
        }
      } else if (ncdims == 2) {
        DW = dimlens[1];
        DH = dimlens[0];
      } else fatal("Unsupported NetCDF variable dimensionality");
      
      // determine DT
      if (DT == 0) DT = input_filenames.size();
      else DT = std::min(DT, input_filenames.size());

      NC_SAFE_CALL( nc_close(ncid) );
#else
      fatal("FTK not compiled with NetCDF.");
#endif
    }
  } 
 
  fprintf(stderr, "SUMMARY\n=============\n");
  fprintf(stderr, "input_filename_pattern=%s\n", input_filename_pattern.c_str());
  fprintf(stderr, "input_format=%s\n", input_format.c_str());
  fprintf(stderr, "output_filename_pattern=%s\n", output_filename_pattern.c_str());
  fprintf(stderr, "output_format=%s\n", output_format.c_str());
  fprintf(stderr, "nd=%d\n", nd);
  fprintf(stderr, "ncdims=%d\n", ncdims);
  fprintf(stderr, "input_variable_name=%s\n", input_variable_name.c_str());
  fprintf(stderr, "DW=%zu\n", DW);
  fprintf(stderr, "DH=%zu\n", DH);
  fprintf(stderr, "DD=%zu\n", DD);
  fprintf(stderr, "DT=%zu\n", DT);
  fprintf(stderr, "threshold=%f\n", threshold);
  fprintf(stderr, "nthreads=%d\n", nthreads);
  fprintf(stderr, "=============\n");

  assert(nd == 2 || nd == 3);
  assert(DT > 0);

  return 0;
}

void track_levelset()
{
  auto *tracker = new ftk::levelset_tracker<>; // ftk::connected_component_tracker<>;
  tracker->set_threshold( threshold );

  for (int current_timestep = 0; current_timestep < DT; current_timestep ++) {
    fprintf(stderr, "current_timestep=%d\n", current_timestep);
    ftk::ndarray<double> field_data = request_timestep(current_timestep);
    // ftk::ndarray<size_t> label_data = threshold_filter<size_t>(field_data);
    // size_t nc = ftk::hoshen_kopelman_2d(label_data);
    // tracker->push_labeled_data_snapshot(label_data.std_vector());
    
    tracker->push_scalar_field_data_snapshot(field_data);
    tracker->advance_timestep();
  }

  tracker->finalize();

  const auto &tg = tracker->get_tracking_graph();
  tg.generate_dot_file("dot");

  delete tracker;
}

int main(int argc, char **argv)
{
  parse_arguments(argc, argv);

  // track_levelset();
  return 0;
}
