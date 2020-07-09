#include <fstream>
#include <mutex>
#include <set>
#include <cassert>
#include "ftk/external/cxxopts.hpp"
#include "ftk/ndarray/synthetic.hh"
#include "ftk/filters/critical_point_tracker_2d_regular.hh"
#include "ftk/filters/critical_point_tracker_3d_regular.hh"
#include "ftk/ndarray.hh"
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
std::string input_filename_pattern, 
  input_format,
  input_dimension,
  input_variable_name, 
  input_variable_name_u, 
  input_variable_name_v, 
  input_variable_name_w;
std::string output_filename,
  output_format;
std::vector<std::string> input_filenames; // assuming each file contains only one timestep, and all files have the exactly same structure
std::string accelerator;
std::string type_filter_str;
size_t DW = 0, DH = 0, DD = 0, DT = 0;
int nthreads = std::thread::hardware_concurrency();
bool verbose = false, demo = false, show_vtk = false, help = false;
bool use_type_filter = false;
unsigned int type_filter = 0;

// determined later
int nd, // dimensionality
    nv; // number of variables; 1 is scalar, otherwise vector
int varid = -1, // only for vti and netcdf
    varid_u = -1,
    varid_v = -1,
    varid_w = -1;
int ncdims = 0;  // Only for netcdf.  
                 // ndndims equals to either nd or (nd+1). 
                 // In the latter case, one of the netcdf dimension 
                 // is time.
int dimids[4] = {-1}; // Only for netcdf
size_t dimlens[4] = {0}; // Only for netcdf

// tracker
ftk::critical_point_tracker_regular* tracker = NULL;

// constants
static const std::set<std::string> set_valid_output_format({str_auto, str_text, str_vtp});


///////////////////////////////
ftk::ndarray<double> request_timestep(int k) // requesting k-th timestep
{
  std::vector<size_t> shape;
  if (nd == 2) {
    if (nv == 1) shape = std::vector<size_t>({DW, DH});
    else shape = std::vector<size_t>({size_t(nv), DW, DH});
  } else {
    if (nv == 1) shape = std::vector<size_t>({DW, DH, DD});
    else shape = std::vector<size_t>({size_t(nv), DW, DH, DD});
  }

  if (demo) {
    if (nd == 2) {
      const double t = DT == 1 ? 0.0 : double(k)/(DT-1);
      return ftk::synthetic_woven_2D<double>(DW, DH, t);
    } else { // nd == 3
      fprintf(stderr, "3D demo case not available.\n");
      assert(false); // TODO: create a 3D demo case
      return ftk::ndarray<double>();
    } 
  } else {
    const std::string filename = input_filenames[k];

    if (input_format == str_float32) {
      ftk::ndarray<float> array32(shape);
      array32.from_binary_file(filename);
      
      ftk::ndarray<double> array(shape);
      array.from_array(array32);

      return array;
    } else if (input_format == str_float64) {
      ftk::ndarray<double> array(shape);
      array.from_binary_file(filename);
      return array;
    } else if (input_format == str_vti) {
      ftk::ndarray<double> array;

      if (input_variable_name.size() > 0) { // all data in one single variable; channels are automatically handled in ndarray
        array.from_vtk_image_data_file(filename, input_variable_name);
      } else { // u, v, w in separate variables
        ftk::ndarray<double> u, v, w;
        u.from_vtk_image_data_file(filename, input_variable_name_u);
        v.from_vtk_image_data_file(filename, input_variable_name_v);
        if (nv > 2)
          w.from_vtk_image_data_file(filename, input_variable_name_w);

        array.reshape(shape);
        for (auto i = 0; i < u.nelem(); i ++) {
          array[i*nv] = u[i];
          array[i*nv+1] = v[i];
          if (nv > 2) array[i*nv+2] = w[i];
        }
      }

      return array;
    } else if (input_format == str_netcdf) {
      ftk::ndarray<double> array;

      if (input_variable_name.size() > 0) { // all data in one single variable; channels are automatically handled in ndarray
        array.from_netcdf(filename, input_variable_name);
      } else { // u, v, w in separate variables
        ftk::ndarray<double> u, v, w;
        u.from_netcdf(filename, input_variable_name_u);
        v.from_netcdf(filename, input_variable_name_v);
        if (nv > 2)
          w.from_netcdf(filename, input_variable_name_w);

        array.reshape(shape);
        for (auto i = 0; i < u.nelem(); i ++) {
          array[i*nv] = u[i];
          array[i*nv+1] = v[i];
          if (nv > 2) array[i*nv+2] = w[i];
        }
      }

      array.reshape(shape); // ncdims may not be equal to nd
      return array;
    } else if (input_format == str_hdf5) {
      // TODO
      assert(false);
      return ftk::ndarray<double>();
    } else {
      assert(false);
      return ftk::ndarray<double>();
    }
  }
}

int parse_arguments(int argc, char **argv)
{
  cxxopts::Options options(argv[0]);
  options.add_options()
    ("i,input", "Input file name pattern: a single file or a series of file, e.g. 'scalar.raw', 'cm1out_000*.nc'", 
     cxxopts::value<std::string>(input_filename_pattern))
    ("demo", "Use synthetic data for demo", cxxopts::value<bool>(demo))
    ("f,input-format", "Input file format (auto|float32|float64|nc|h5|vti)", 
     cxxopts::value<std::string>(input_format)->default_value(str_auto))
    ("dim", "Spatial dimensionality of data (auto|2|3)", 
     cxxopts::value<std::string>(input_dimension)->default_value(str_auto))
    ("nvar", "Number of variables",
     cxxopts::value<int>(nv)->default_value("0"))
    ("w,width", "Width", cxxopts::value<size_t>(DW))
    ("h,height", "Height", cxxopts::value<size_t>(DH))
    ("d,depth", "Depth", cxxopts::value<size_t>(DH))
    ("n,timesteps", "Number of timesteps", cxxopts::value<size_t>(DT))
    ("var", "Variable name (only for NetCDF, HDF5, and VTK)", 
     cxxopts::value<std::string>(input_variable_name))
    ("var-u", "Variable name for u-component",
     cxxopts::value<std::string>(input_variable_name_u))
    ("var-v", "Variable name for v-component",
     cxxopts::value<std::string>(input_variable_name_v))
    ("var-w", "Variable name for w-component",
     cxxopts::value<std::string>(input_variable_name_w))
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

  if (input_variable_name.size() != 0 &&
      input_variable_name_u.size() + 
      input_variable_name_v.size() + 
      input_variable_name_w.size() != 0) {
    fprintf(stderr, "input_variable_name=%s\n", input_variable_name.c_str());
    fprintf(stderr, "input_variable_name_u=%s\n", input_variable_name_u.c_str());
    fprintf(stderr, "input_variable_name_v=%s\n", input_variable_name_v.c_str());
    fprintf(stderr, "input_variable_name_w=%s\n", input_variable_name_w.c_str());
    fatal("Cannot specify both `--var' and `--var-w|--var-v|--var-w' simultanuously");
  }

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

  // processing input
  if (demo) {
    if (input_filename_pattern.size())
      fprintf(stderr, "WARN: '--input' ignored.\n");
    if (input_format != "auto")
      fprintf(stderr, "WARN: '--input-format' ignored.\n");

    if (input_dimension == "3") nd = 3;
    else nd = 2;

    if (nv == 0) nv = 1; // use scalar data for demo
    
    if (input_variable_name.size())
      warn("--var ignored.");
    if (input_variable_name_u.size())
      warn("--var-u ignored.");
    if (input_variable_name_v.size())
      warn("--var-v ignored.");
    if (input_variable_name_w.size())
      warn("--var-w ignored.");

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
      if (input_variable_name.size() + input_variable_name_u.size() 
          + input_variable_name_v.size() + input_variable_name_w.size())
        warn("Ignoring variable names.");

      // determine nv
      if (nv == 0) // auto
        fatal("Unable to determine the number of variables"); // TOOD: determine nv by file size
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

      // determine variable names and array indices
      if (input_variable_name.size() == 0) {
        if (input_variable_name_u.size() == 0) { // data from single unamed variable
          varid = 0;
          input_variable_name = image->GetPointData()->GetArrayName(varid);
          // fprintf(stderr, "Using input variable: %s\n", input_variable_name.c_str());
        } else { // data from multiple named variables
          // nv = nd;
          image->GetPointData()->GetArray(input_variable_name_u.c_str(), varid_u);
          if (varid_u < 0) fatal("Cannot find variable name for u.");
          image->GetPointData()->GetArray(input_variable_name_v.c_str(), varid_v);
          if (varid_v < 0) fatal("Cannot find variable name for v.");
          if (nd == 3) {
            image->GetPointData()->GetArray(input_variable_name_w.c_str(), varid_w);
            if (varid_w < 0) fatal("Cannot find variable name for w.");
          }
        }
      } else { // data from single named variable
        image->GetPointData()->GetArray(input_variable_name.c_str(), varid);
      }

      // determine nv
      if (nv == 0) {
        if (varid >= 0) 
          nv = image->GetPointData()->GetArray(varid)->GetNumberOfComponents();
        else 
          nv = nd; // TODO: make sure each of u,v,w has single component.
        // TODO: sanity check
      } else 
        fatal("Cannot override number of variables in VTI.");

      // determine DT
      if (DT == 0) DT = input_filenames.size();
      else DT = std::min(DT, input_filenames.size());
#else
      fatal("FTK not compiled with VTK.");
#endif
    } else if (input_format == str_netcdf) {
#if FTK_HAVE_NETCDF
      if (input_variable_name.size() +
          input_variable_name_u.size() + 
          input_variable_name_v.size() + 
          input_variable_name_w.size() == 0) 
        fatal("Variable name missing for NetCDF files.");

      int ncid, my_varid;
      NC_SAFE_CALL( nc_open(input_filenames[0].c_str(), NC_NOWRITE, &ncid) );

      if (input_variable_name.size() > 0) { // single variable
        NC_SAFE_CALL( nc_inq_varid(ncid, input_variable_name.c_str(), &varid) );
        my_varid = varid;
        nv = 1; // TODO: netcdf multicomponent variables
      } else {
        NC_SAFE_CALL( nc_inq_varid(ncid, input_variable_name_u.c_str(), &varid_u) );
        NC_SAFE_CALL( nc_inq_varid(ncid, input_variable_name_v.c_str(), &varid_v) );
        if (input_variable_name_w.size() > 0) {
          NC_SAFE_CALL( nc_inq_varid(ncid, input_variable_name_w.c_str(), &varid_w) );
          nv = 3;
        } else
          nv = 2;
        my_varid = varid_u;
        // TODO: check if u, v, w have the same dimensions
      }
      NC_SAFE_CALL( nc_inq_varndims(ncid, my_varid, &ncdims) );
      
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
      NC_SAFE_CALL( nc_inq_vardimid(ncid, my_varid, dimids) );
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
  fprintf(stderr, "output_filename=%s\n", output_filename.c_str());
  fprintf(stderr, "output_format=%s\n", output_format.c_str());
  fprintf(stderr, "nd=%d\n", nd);
  fprintf(stderr, "ncdims=%d\n", ncdims);
  fprintf(stderr, "nv=%d\n", nv);
  if (input_variable_name.size())
    fprintf(stderr, "input_variable_name=%s\n", input_variable_name.c_str());
  if (input_variable_name_u.size())
    fprintf(stderr, "input_variable_name_u=%s\n", input_variable_name_u.c_str());
  if (input_variable_name_v.size())
    fprintf(stderr, "input_variable_name_v=%s\n", input_variable_name_v.c_str());
  if (input_variable_name_w.size())
    fprintf(stderr, "input_variable_name_w=%s\n", input_variable_name_w.c_str());
  fprintf(stderr, "DW=%zu\n", DW);
  fprintf(stderr, "DH=%zu\n", DH);
  fprintf(stderr, "DD=%zu\n", DD);
  fprintf(stderr, "DT=%zu\n", DT);
  fprintf(stderr, "type_filter=%s\n", type_filter_str.c_str());
  fprintf(stderr, "nthreads=%d\n", nthreads);
  fprintf(stderr, "=============\n");

  assert(nd == 2 || nd == 3);
  assert(nv == 1 || nv == 2 || nv == 3);
  assert(DT > 0);

  return 0;
}

void track_critical_points()
{
  if (nd == 2) {
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

  int current_timestep = 0;
  while (1) {
    ftk::ndarray<double> field_data = request_timestep(current_timestep);
    if (nv == 1) // scalar field
      tracker->push_scalar_field_snapshot(field_data);
    else // vector field
      tracker->push_vector_field_snapshot(field_data);
     
    if (current_timestep == DT - 1) {
      tracker->update_timestep();
      break;
    }
    else if (current_timestep != 0) // need to push two timestep before one can advance timestep
      tracker->advance_timestep();
    current_timestep ++;
  }

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
