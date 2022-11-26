#ifndef _FTK_REGULAR_ARRAY_STREAM_HH
#define _FTK_REGULAR_ARRAY_STREAM_HH

#include <ftk/object.hh>
#include <fstream>
#include <chrono>
#include <condition_variable>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/filters/streaming_filter.hh>
#include <ftk/external/json.hh>
#include <ftk/utils/scatter.hh>
#include <ftk/utils/bcast.hh>

#if FTK_HAVE_VTK
#include <vtkResampleToImage.h>
#endif

namespace ftk {
using nlohmann::json;

template <typename T=double>
struct ndarray_stream : public object {
  ndarray_stream(diy::mpi::communicator comm = MPI_COMM_WORLD);
  ndarray_stream(const std::string adios2_config_filename, //  = "adios2.xml",
      const std::string adios2_io_name, //  = "SimulationOutput",
      diy::mpi::communicator comm); //  = MPI_COMM_WORLD);
  ~ndarray_stream();

  void configure(const json& j) {set_input_source_json(j);}
  // JSON specifications:
  // required fields: 
  //  - type, string.  Must be one of "synthetic" and "file".  This field may be ommited
  //    if `format' is given.
  //  - name (required if type is synthetic), string.  Must be one of the follows: "woven", 
  //    "double_gyre", "merger_2d".
  // optional fields:
  //  - filenames (required if type is file), string.  The list of filenames will be 
  //    determined by glob(3)
  //  - mesh_filename (required for unstructured mesh, must be in vtu format in the 
  //    current version), string.
  //  - format (required if type is file and format is float32/float64), string.  If not 
  //    given, the format will be determined by the filename extension.  The value of this 
  //    field must be one of the follows: vti, vtu, nc, h5, float32, float64.
  //  - variables (required if format is nc/h5, optional for vti/vtu), array of strings.
  //    - the number of components is the length of the array.
  //    - if not given, the defaulat value is ["scalar"]
  //  - components (to be determined), array of number of components per variable
  //  - dimensions (required if format is floaot32/float64), array of integers, e.g.  
  //    [width, height, depth]
  //  - n_timesteps, integer.  The default is 32 for synthetic data; the number can be 
  //    automatically derived from file; the number can be automatically derived from files
  //  - perturbation, number.  Add gaussian perturbation to the data
  //  - clamp, array of two numbers (min, max).  Clamp the range of the input data 
  //    with the given min and max values

  void set_input_source_json_file(const std::string& filename);
  void set_input_source_json(const json& j_);
  const json& get_json() const {return j;}
  
  void start();
  void finish();

  void set_callback(std::function<void(int, const ndarray<T>&)> f) {callback = f;}

  size_t n_variables() const { return j["variables"].size(); }
  size_t n_components() const {
    size_t n = 0;
    for (int i = 0; i < j["components"].size(); i ++)
      n += j["components"][i].template get<size_t>();
    return n;
  }
  size_t n_dimensions() const {
    // if (j.contains("nc_has_unlimited_time_dimension")) return j["dimensions"].size() - 1;
    // else return j["dimensions"].size(); 
    return j["dimensions"].size(); 
  }
  size_t n_timesteps() const {
    if (j.contains("n_timesteps")) return j["n_timesteps"];
    else return std::numeric_limits<size_t>::max();
  }

  std::vector<size_t> shape() const;

protected:
  ndarray<T> request_timestep_file(int k);
  ndarray<T> request_timestep_file_nc(int k);
  ndarray<T> request_timestep_file_pnc(int k);
  ndarray<T> request_timestep_file_vti(int k);
  ndarray<T> request_timestep_file_vtu(int k);
  ndarray<T> request_timestep_file_h5(int k);
  ndarray<T> request_timestep_file_bp3(int k);
  ndarray<T> request_timestep_file_bp4(int k);
  template <typename T1> ndarray<T> request_timestep_file_binary(int k);

  ndarray<T> request_timestep_synthetic(int k);
  ndarray<T> request_timestep_synthetic_woven(int k);
  ndarray<T> request_timestep_synthetic_moving_extremum_2d(int k);
  ndarray<T> request_timestep_synthetic_moving_extremum_3d(int k);
  ndarray<T> request_timestep_synthetic_moving_ramp_3d(int k);
  ndarray<T> request_timestep_synthetic_moving_dual_ramp_3d(int k);
  ndarray<T> request_timestep_synthetic_double_gyre(int k);
  ndarray<T> request_timestep_synthetic_merger_2d(int k);
  ndarray<T> request_timestep_synthetic_volcano_2d(int k);
  ndarray<T> request_timestep_synthetic_tornado(int k);

  void modified_callback(int, const ndarray<T>&);

public:
  bool is_partial_read_supported() const;

  void set_part(const lattice& e) { 
    part = true; ext = e; 
    std::cerr << "partial domain: " << e << std::endl; 
  }

protected:
  bool part = false;
  lattice ext;

protected:
  json j; // configs, metadata, and everything

  std::function<void(int, const ndarray<T>&)> callback;

  streaming_filter<ndarray<T>, T> temporal_filter;

protected: // adios2
#if FTK_HAVE_ADIOS2
  adios2::ADIOS adios;
  adios2::IO adios_io;
  adios2::Engine adios_reader;
  std::vector<adios2::Variable<T>> adios_vars;
  int current_file_id = 0;
#endif
  std::string adios2_config_filename, 
              adios2_io_name = "BPReader";

private:
  static bool ends_with(std::string const & value, std::string const & ending)
  {
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
  }
};

template <typename T>
ndarray_stream<T>::ndarray_stream(diy::mpi::communicator comm)
#if ADIOS2_USE_MPI
  : adios(comm)
#endif
{
#if FTK_HAVE_ADIOS2
  adios_io = adios.DeclareIO(adios2_io_name);
#endif
}

template <typename T>
ndarray_stream<T>::ndarray_stream(
    const std::string adios2_config_filename, 
    const std::string adios2_io_name,
    diy::mpi::communicator comm)
// #if FTK_HAVE_ADIOS2
#if ADIOS2_USE_MPI
  : adios(adios2_config_filename, comm, adios2::DebugON)
#endif
{
  fprintf(stderr, "creating stream...\n");
#if FTK_HAVE_ADIOS2
  this->adios2_config_filename = adios2_config_filename;
  this->adios2_io_name = adios2_io_name;
  adios_io = adios.DeclareIO(adios2_io_name);
#endif
}

template <typename T>
ndarray_stream<T>::~ndarray_stream()
{
}

template <typename T>
void ndarray_stream<T>::set_input_source_json_file(const std::string& filename)
{
  std::ifstream t(filename);
  std::string str((std::istreambuf_iterator<char>(t)),
                   std::istreambuf_iterator<char>());
  json j = json::parse(str);
  set_input_source_json(j);
  t.close();
}

template <typename T>
void ndarray_stream<T>::set_input_source_json(const json& j_)
{
  j = j_;
  // std::cerr << j << std::endl;
 
  j["adios2_config"] = adios2_config_filename;
  j["adios2_name"] = adios2_io_name;

  if (!j.contains("type") && j.contains("format"))
    j["type"] = "file";

  // check if missing dimensions
  bool missing_dimensions = false;
  if (j.contains("dimensions")) {
    if (j["dimensions"].is_array()) {
      for (int i = 0; i < j["dimensions"].size(); i ++) {
        if (j["dimensions"][i].is_number()) {
          // OK
        } else
          fatal("invalid array dimensions (1)");
      }
    } else 
      fatal("invalid array dimensions (2)");
  } else 
    missing_dimensions = true;

  // check bounds; 
  if (j.contains("bounds")) { // TODO: check if inputs are an array or a string
    auto jb = j["bounds"];
    if (jb.is_string()) { // convert a comma-separeted string to an array
      std::vector<double> mybounds;
      if (jb.is_string()) { // parse comma separated bounds
        auto strs = ftk::split(jb.template get<std::string>(), ",");
        for (const auto str : strs) {
          mybounds.push_back(std::stod(str));
          // fprintf(stderr, "%f\n", std::stod(str));
        }
      }
      j["bounds"] = mybounds;
    } else if (jb.is_array()) { // check if each element is a number; the maximum length of the array should be 6
      if (jb.size() > 6)
        warn("bounds cannot exceed 6 numbers");

      for (int i = 0; i < jb.size(); i ++)
        if (!jb[i].is_number())
          fatal("bounds must be numbers");
    } else 
      fatal("bounds must be an array or a comma-seprated string");
  }

  // check if missing variables
  bool missing_variables = false;
  if (j.contains("variables")) { 
    if (j["variables"].is_array()) { 
      for (const auto &v : j["variables"]) {
        if (v.is_string()) {
          // OK
        } else
          fatal("invalid variable name");
      }
    } else 
      fatal("invalid variable list");
  } else missing_variables = true;

  if (j.contains("type")) {
    if (j["type"] == "synthetic") {
      int default_nd;
      int default_dims[3] = {32, 32, 32};
      int default_n_timesteps = 32;
      std::vector<int> default_components = {1};

      if (j.contains("name")) {
        if (j["name"] == "woven") {
          if (missing_variables) 
            j["variables"] = {"scalar"};
          default_nd = 2;
          if (!j.contains("x0")) j["scalaring_factor"] = 15.0;
        } else if (j["name"] == "moving_extremum_2d") {
          if (missing_variables) 
            j["variables"] = {"scalar"};
          default_nd = 2;
          default_dims[0] = 21; 
          default_dims[1] = 21;
          default_dims[2] = 21;

          if (j.contains("x0")) {
            if (j["x0"].is_array()) {
              if (j["x0"].size() == 2) {
                for (int i = 0; i < 2; i ++)
                  if (j["x0"][i].is_number()) {
                    // OK
                  } else 
                    fatal("invalid x0 coordinates");
              } else
                fatal("invalid x0");
            } else 
              fatal("invalid x0");
          } else 
            j["x0"] = {10.0, 10.0};

          if (j.contains("dir")) {
            if (j["dir"].is_array()) {
              if (j["dir"].size() == 2) {
                for (int i = 0; i < 2; i ++)
                  if (j["dir"][i].is_number()) {
                    // OK
                  } else 
                    fatal("invalid dir");
              } else
                fatal("invalid dir");
            } else 
              fatal("invalid dir");
          } else 
            // j["dir"] = {0.5, 0.5};
            j["dir"] = {0.1, 0.1};
        } else if (j["name"] == "moving_extremum_3d") {
          if (missing_variables) 
            j["variables"] = {"scalar"};
          default_nd = 3;
          default_dims[0] = 21; 
          default_dims[1] = 21;
          default_dims[2] = 21;

          if (j.contains("x0")) {
            if (j["x0"].is_array()) {
              if (j["x0"].size() == 3) {
                for (int i = 0; i < 3; i ++)
                  if (j["x0"][i].is_number()) {
                    // OK
                  } else 
                    fatal("invalid x0 coordinates");
              } else
                fatal("invalid x0");
            } else 
              fatal("invalid x0");
          } else 
            j["x0"] = {10, 10, 10};

          if (j.contains("dir")) {
            if (j["dir"].is_array()) {
              if (j["dir"].size() == 3) {
                for (int i = 0; i < 3; i ++)
                  if (j["dir"][i].is_number()) {
                    // OK
                  } else 
                    fatal("invalid dir");
              } else
                fatal("invalid dir");
            } else 
              fatal("invalid dir");
          } else 
            j["dir"] = {0.1, 0.11, 0.1};
          
        } else if (j["name"] == "moving_ramp_3d") {
          if (missing_variables) 
            j["variables"] = {"scalar"};
          default_nd = 3;
          default_dims[0] = 21; 
          default_dims[1] = 21;
          default_dims[2] = 21;
          
          if (j.contains("x0")) {
            if (j["x0"].is_number()) { // OK
            } else fatal("invalid x0");
          } else j["x0"] = 10;

          if (j.contains("rate")) {
            if (j["rate"].is_number()) { // OK
            } else fatal("invalid rate");
          } else j["rate"] = 0.1;
        } else if (j["name"] == "moving_dual_ramp_3d") {
          if (missing_variables) 
            j["variables"] = {"scalar"};
          default_nd = 3;
          default_dims[0] = 21; 
          default_dims[1] = 21;
          default_dims[2] = 21;
          
          if (j.contains("x0")) {
            if (j["x0"].is_number()) { // OK
            } else fatal("invalid x0");
          } else j["x0"] = 10;

          if (j.contains("rate")) {
            if (j["rate"].is_number()) { // OK
            } else fatal("invalid rate");
          } else j["rate"] = 0.7;

          if (j.contains("offset")) {
            if (j["offset"].is_number()) { // OK
            } else fatal("invalid offset");
          } else j["offset"] = 2.0;
        } else if (j["name"] == "volcano_2d") {
          if (missing_variables) 
            j["variables"] = {"scalar"};
          default_nd = 2;
          default_dims[0] = 21; 
          default_dims[1] = 21;
          default_dims[2] = 21;
          
          if (j.contains("x0")) {
            if (j["x0"].is_array()) {
              if (j["x0"].size() == 2) {
                for (int i = 0; i < 2; i ++)
                  if (j["x0"][i].is_number()) {
                    // OK
                  } else 
                    fatal("invalid x0 coordinates");
              } else
                fatal("invalid x0");
            } else 
              fatal("invalid x0");
          } else 
            j["x0"] = {10.0, 10.0};
          
          if (j.contains("radius")) {
            if (j["radius"].is_number()) { // OK
            } else fatal("invalid radius");
          } else 
            j["radius"] = 3.0;
          
        } else if (j["name"] == "double_gyre") {
          default_nd = 2;
          default_dims[0] = 64; 
          default_dims[1] = 32;
          default_n_timesteps = 50;
          j["variables"] = {"u", "v", "w"};
          j["components"] = {1, 1, 1};
          if (j.contains("time_scale")) {
            if (j["time_scale"].is_number()) { // ok
            } else fatal("invalid time_scale");
          } else j["time_scale"] = 0.1;
        } else if (j["name"] == "merger_2d") {
          j["variables"] = {"scalar"};
          default_nd = 2;
          default_dims[0] = 32; 
          default_dims[1] = 32;
          default_n_timesteps = 100;
        } else if (j["name"] == "tornado") {
          default_nd = 3;
          j["variables"] = {"u", "v", "w"};
          j["components"] = {1, 1, 1};
        } else {
          std::cerr << "synthetic case name: " << j["name"] << std::endl;
          fatal("synthetic case not available.");
        }
      } else fatal("synthetic case name not given.");
     
      if (missing_dimensions) {
        std::vector<int> dims;
        for (int i = 0; i < default_nd; i ++)
          dims.push_back(default_dims[i]);
        j["dimensions"] = dims;
      }

      if (!j.contains("components"))
        j["components"] = default_components;

      // if (j.contains("n_timesteps")) assert(j["n_timesteps"] != 0);
      if (!j.contains("n_timesteps"))
        j["n_timesteps"] = default_n_timesteps;
    } else if (j["type"] == "file") {
      if (j.contains("filenames")) {
        if (j["filenames"].is_array()) {
          // j["n_timesteps"] = j["filenames"].size(); // TODO: we are assuming #timesteps = #filenames
          j["n_timesteps"] = j["filenames"].size(); // TODO: we are assuming #timesteps = #filenames
        } else {
          auto filenames = glob(j["filenames"]);
          if (filenames.empty()) {
            if (ends_with(j["filenames"], "bp")) {
              // OK, input are adios streams
              filenames.push_back( j["filenames"] ); // one single "filename"
              j["format"] = "bp4"; // the input format must be adios2/bp4
              
              const std::string filename0 = j["filenames"];
              j["filenames"] = { filename0 };
            } else {
              fatal("unable to find matching filename(s).");
            }
          } else 
            j["filenames"] = filenames;
         
#if 0
          if (j.contains("n_timesteps") && j["n_timesteps"].is_number()) { 
            j["n_timesteps"] = std::min(j["n_timesteps"].template get<size_t>(), j["filenames"].size());
            filenames.resize(j["n_timesteps"]);
          }
          else 
            j["n_timesteps"] = j["filenames"].size();
#endif
        }
        // std::cerr << j << std::endl;
        const std::string filename0 = j["filenames"][0];

        if (!j.contains("format")) { // probing file format
          const auto ext = file_extension(filename0);
          
          if (ext == FILE_EXT_VTI) j["format"] = "vti";
          else if (ext == FILE_EXT_VTU) j["format"] = "vtu";
          else if (ext == FILE_EXT_PVTU) j["format"] = "pvtu";
          else if (ext == FILE_EXT_NETCDF) j["format"] = "nc";
          else if (ext == FILE_EXT_HDF5) j["format"] = "h5";
          else if (ext == FILE_EXT_BP) { // need to further distinguish if input is bp3 or bp4
            if (is_directory(filename0)) j["format"] = "bp4";
            else j["format"] = "bp3";
          }
          else fatal(FTK_ERR_FILE_UNRECOGNIZED_EXTENSION);
        }

        if (j["format"] == "float32" || j["format"] == "float64") {
          if (missing_dimensions) 
            fatal("missing dimensions.");
        
          if (missing_variables) { // add one single component in the array
            const std::string var("scalar");
            const std::vector<std::string> vars = {var};
            j["variables"] = vars;
            j["components"] = {1};
          } else {
            const size_t nv = j["variables"].size();
            const std::vector<int> ones(nv, 1);
            j["components"] = ones;
          }
          
          if (j.contains("n_timesteps") && j["n_timesteps"].is_number())
            j["n_timesteps"] = std::min(j["n_timesteps"].template get<size_t>(), j["filenames"].size());
          else 
            j["n_timesteps"] = j["filenames"].size();
        } else if (j["format"] == "vti") {
#if FTK_HAVE_VTK
          if (is_root_proc()) {
            vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
            reader->SetFileName(filename0.c_str());
            reader->Update();

            vtkSmartPointer<vtkImageData> image = reader->GetOutput();

            if (j.contains("dimensions")) 
              warn("ignoring dimensions");
            int imageNd = image->GetDataDimension();
            if (imageNd == 2)
              j["dimensions"] = {image->GetDimensions()[0], image->GetDimensions()[1]};
            else 
              j["dimensions"] = {image->GetDimensions()[0], image->GetDimensions()[1], image->GetDimensions()[2]};
         
            if (j.contains("bounds"))
              warn("input bounds will be overridden");

            std::vector<double> bounds(6, 0);
            image->GetBounds(&bounds[0]);
            j["bounds"] = bounds;

            if (missing_variables) {
              const std::string var = image->GetPointData()->GetArrayName(0);
              j["variables"] = {var};
            }

            // determine number of components per var
            std::vector<int> components;
            for (int i = 0; i < j["variables"].size(); i ++) {
              const std::string var = j["variables"][i];
              vtkSmartPointer<vtkDataArray> da = image->GetPointData()->GetArray( var.c_str() );
              if (!da) fatal(FTK_ERR_VTK_VARIABLE_NOT_FOUND);
              const int nc = da->GetNumberOfComponents();
              components.push_back(nc);
            }
            j["components"] = components;
            
            // determine number timesteps
            if (j.contains("n_timesteps") && j["n_timesteps"].is_number())
              j["n_timesteps"] = std::min(j["n_timesteps"].template get<size_t>(), j["filenames"].size());
            else 
              j["n_timesteps"] = j["filenames"].size();
          }
          
          diy::mpi::bcastj(comm, j, get_root_proc());
#else
          fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
        } else if (j["format"] == "vtu" || j["format"] == "pvtu" || j["format"] == "vtu_resample" || j["format"] == "pvtu_resample") {
#if FTK_HAVE_VTK
          if (is_root_proc()) {
            bool resample = false;
            if (j["format"] == "vtu_resample" || j["format"] == "pvtu_resample") 
              resample = true;

            vtkSmartPointer<vtkUnstructuredGrid> grid;
            if (j["format" ] == "vtu" || j["format"] == "vtu_resample") {
              vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
              reader->SetFileName(filename0.c_str());
              reader->Update();
              grid = reader->GetOutput();
            } else if (j["format"] == "pvtu" || j["format"] == "pvtu_resample") {
              vtkSmartPointer<vtkXMLPUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();
              reader->SetFileName(filename0.c_str());
              reader->Update();
              grid = reader->GetOutput();
            }
           
            if (resample) {
              if (missing_dimensions)
                fatal("missing dimensions.");

              // if resample bounds are missing, will use auto bounds
            } else {
              if (j.contains("dimensions")) 
                warn("ignoring dimensions");
              j["dimensions"] = {0}; // workaround
            }
            
            if (missing_variables) {
              const std::string var = grid->GetPointData()->GetArrayName(0);
              j["variables"] = {var};
            }
            
            // determine number of components per var
            std::vector<int> components;
            for (int i = 0; i < j["variables"].size(); i ++) {
              const std::string var = j["variables"][i];
              vtkSmartPointer<vtkDataArray> da = grid->GetPointData()->GetArray( var.c_str() );
              if (!da) fatal(FTK_ERR_VTK_VARIABLE_NOT_FOUND);
              const int nc = da->GetNumberOfComponents();
              components.push_back(nc);
            }
            j["components"] = components;
            
            // determine number timesteps
            if (j.contains("n_timesteps") && j["n_timesteps"].is_number())
              j["n_timesteps"] = std::min(j["n_timesteps"].template get<size_t>(), j["filenames"].size());
            else 
              j["n_timesteps"] = j["filenames"].size();
          }

          diy::mpi::bcastj(comm, j, get_root_proc());
#else 
          fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
        } else if (j["format"] == "nc") {
#if FTK_HAVE_NETCDF
          if (missing_variables) 
            fatal(FTK_ERR_NETCDF_MISSING_VARIABLE);

          const int nv = j["variables"].size();
          int ncid, ncdims, nd, varids[nv];

#if NC_HAS_PARALLEL
          int rtn = nc_open_par(filename0.c_str(), NC_NOWRITE, comm, MPI_INFO_NULL, &ncid);
          if (rtn != NC_NOERR)
            NC_SAFE_CALL( nc_open(filename0.c_str(), NC_NOWRITE, &ncid) );
#else
          NC_SAFE_CALL( nc_open(filename0.c_str(), NC_NOWRITE, &ncid) );
#endif

          for (int i = 0; i < nv; i ++) {
            const std::string var = j["variables"][i];
            NC_SAFE_CALL( nc_inq_varid(ncid, var.c_str(), &varids[i]) );
          }
          NC_SAFE_CALL( nc_inq_varndims(ncid, varids[0], &ncdims) ); // assuming all variables have the same dimensions

          // components
          const std::vector<int> components(nv, 1);
          j["components"] = components;

          // determin spatial dimensions
          int dimids[4]; 
          size_t dimlens[4];
          NC_SAFE_CALL( nc_inq_vardimid(ncid, varids[0], dimids) );
          for (int i = 0; i < ncdims; i ++) 
            NC_SAFE_CALL( nc_inq_dimlen(ncid, dimids[i], &dimlens[i]) );

          int nt;
          int unlimited_recid;
          NC_SAFE_CALL( nc_inq_unlimdim(ncid, &unlimited_recid) );
          NC_SAFE_CALL( nc_close(ncid) );
          if (unlimited_recid >= 0) {
            j["nc_has_unlimited_time_dimension"] = true;

            // check timesteps per file
            const int nf = j["filenames"].size();
            std::vector<int> timesteps_per_file(nf), first_timestep_per_file(nf);
            int total_timesteps = 0;
            for (int i = 0; i < nf; i ++) {
              const std::string filename = j["filenames"][i];
              size_t nt;
              NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
              NC_SAFE_CALL( nc_inq_dimlen(ncid, dimids[0], &nt) );
              NC_SAFE_CALL( nc_close(ncid) );
              timesteps_per_file[i] = nt;
              first_timestep_per_file[i] = total_timesteps;
              total_timesteps += nt;
            }
            j["timesteps_per_file"] = timesteps_per_file;
            j["first_timestep_per_file"] = first_timestep_per_file;
          
            nt = std::accumulate(timesteps_per_file.begin(), timesteps_per_file.end(), 0);
          } else {
            nt = j["filenames"].size();
          }
            
          if (j.contains("n_timesteps") && j["n_timesteps"].is_number())
            j["n_timesteps"] = std::min(j["n_timesteps"].template get<int>(), nt);
          else 
            j["n_timesteps"] = nt;
        
          std::vector<int> given_dimensions;
          if (j.contains("dimensions")) {
            auto dims = j["dimensions"];
            for (int i = 0; i < dims.size(); i ++)
              given_dimensions.push_back(dims[i]);
            // fprintf(stderr, "#given_dimensions=%zu\n", given_dimensions.size());
            // warn("ignorning dimensions");
          }

          if (given_dimensions.size() == 0) {
            if (ncdims == 4) { // 3 spatial dims + 1 time dimension; nd not required
              j["dimensions"] = {dimlens[3], dimlens[2], dimlens[1]};
            } else if (ncdims == 3) {
              if (unlimited_recid >= 0)
                j["dimensions"] = {dimlens[2], dimlens[1]};
              else 
                j["dimensions"] = {dimlens[2], dimlens[1], dimlens[0]};
            } else if (ncdims == 2) {
              if (unlimited_recid >= 0)
                j["dimensions"] = {dimlens[1]};
              else 
                j["dimensions"] = {dimlens[1], dimlens[0]};
            } else 
              fatal("unsupported netcdf variable dimensionality");
          }
#else
          fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
        } else if (j["format"] == "pnc") { // parallel-netcdf
#if FTK_HAVE_PNETCDF
          if (missing_variables)
            fatal(FTK_ERR_NETCDF_MISSING_VARIABLE);

          const int nv = j["variables"].size();
          int ncid, ncdims, nd, varids[nv];

          PNC_SAFE_CALL( ncmpi_open(comm, filename0.c_str(), NC_NOWRITE, MPI_INFO_NULL, &ncid) );
          for (int i = 0; i < nv; i ++) {
            const std::string var = j["variables"][i];
            PNC_SAFE_CALL( ncmpi_inq_varid(ncid, var.c_str(), &varids[i]) );
          }
          PNC_SAFE_CALL( ncmpi_inq_varndims(ncid, varids[0], &ncdims) ); // assuming all variables have the same dimensions
          
          // components
          const std::vector<int> components(nv, 1);
          j["components"] = components;
          
          // determin spatial dimensions
          int dimids[4]; 
          MPI_Offset dimlens[4];
          PNC_SAFE_CALL( ncmpi_inq_vardimid(ncid, varids[0], dimids) );
          for (int i = 0; i < ncdims; i ++) 
            PNC_SAFE_CALL( ncmpi_inq_dimlen(ncid, dimids[i], &dimlens[i]) );

          int nt;
          int unlimited_recid;
          PNC_SAFE_CALL( ncmpi_inq_unlimdim(ncid, &unlimited_recid) );
          PNC_SAFE_CALL( ncmpi_close(ncid) );
          if (unlimited_recid >= 0) {
            j["nc_has_unlimited_time_dimension"] = true;

            // check timesteps per file
            const int nf = j["filenames"].size();
            std::vector<int> timesteps_per_file(nf), first_timestep_per_file(nf);
            int total_timesteps = 0;
            for (int i = 0; i < nf; i ++) {
              const std::string filename = j["filenames"][i];
              size_t nt;
              PNC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
              PNC_SAFE_CALL( nc_inq_dimlen(ncid, dimids[0], &nt) );
              PNC_SAFE_CALL( nc_close(ncid) );
              timesteps_per_file[i] = nt;
              first_timestep_per_file[i] = total_timesteps;
              total_timesteps += nt;
            }
            j["timesteps_per_file"] = timesteps_per_file;
            j["first_timestep_per_file"] = first_timestep_per_file;
          
            nt = std::accumulate(timesteps_per_file.begin(), timesteps_per_file.end(), 0);
          } else {
            nt = j["filenames"].size();
          }

#else
          fatal(FTK_ERR_NOT_BUILT_WITH_PNETCDF);
#endif
        } else if (j["format"] == "h5") {
#if FTK_HAVE_HDF5
          if (missing_variables)
            fatal("missing variables for h5");

          const std::string varname0 = j["variables"][0];
          auto fid = H5Fopen(filename0.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
          auto did = H5Dopen2(fid, varname0.c_str(), H5P_DEFAULT);
          auto sid = H5Dget_space(did);
          const int h5ndims = H5Sget_simple_extent_ndims(sid);
          hsize_t h5dims[h5ndims];
          H5Sget_simple_extent_dims(sid, h5dims, NULL);

          std::vector<size_t> dims(h5ndims);
          for (auto i = 0; i < h5ndims; i ++)
            dims[i] = h5dims[i];
          std::reverse(dims.begin(), dims.end());
          
          if (j.contains("dimensions"))
            warn("ignoring dimensions");
          j["dimensions"] = dims;

          H5Fclose(fid);
          
          // j["n_timesteps"] = j["filenames"].size(); // TODO: we are assuming #timesteps = #filenames
          if (j.contains("n_timesteps") && j["n_timesteps"].is_number())
            j["n_timesteps"] = std::min(j["n_timesteps"].template get<size_t>(), j["filenames"].size());
          else 
            j["n_timesteps"] = j["filenames"].size();
          
          // components
          const size_t nv = j["variables"].size();
          const std::vector<int> components(nv, 1);
          j["components"] = components;
#else
          fatal(FTK_ERR_NOT_BUILT_WITH_HDF5);
#endif
        } else if (j["format"] == "bp4") {
#if FTK_HAVE_ADIOS2
          if (missing_variables)
            fatal("missing variables for bp4");

          // adios2::ADIOS adios(comm);
          // adios2::IO io = adios.DeclareIO("BPReader");
          // adios2::Engine reader = io.Open(filename0, adios2::Mode::Read);
          adios_reader = adios_io.Open(filename0, adios2::Mode::Read);
          while (1) {
            auto status = adios_reader.BeginStep( adios2::StepMode::Read, 10.f );
            if (status == adios2::StepStatus::NotReady) {
              usleep(100000);
              continue;
            } else if (status == adios2::StepStatus::OK) {
              break;
            } else fatal("adios2 error");
          }

          for (int i = 0; i < j["variables"].size(); i ++) {
            const std::string varname = j["variables"][i];
            adios2::Variable<T> var = adios_io.template InquireVariable<T>(varname);
            if (!var) {
              fatal("variable not found in bp4");
            } else {
              adios_vars.push_back( var );
            }
          }

          std::vector<size_t> dims(adios_vars[0].Shape());
          std::reverse(dims.begin(), dims.end()); 
          if (j.contains("dimensions"))
            warn("ignoring bp4 dimensions");
          j["dimensions"] = dims;

          // const size_t nsteps0 = adios_vars[0].Steps();
          // fprintf(stderr, "nsteps0=%zu\n", nsteps0);
          // reader.Close();

          if (j.contains("n_timesteps") && j["n_timesteps"].is_number()) {
            // number of timesteps explicitly specified
            // j["n_timesteps"] = std::min(j["n_timesteps"].template get<size_t>(), j["filenames"].size());
          } else {
            // if (nsteps0 > 1) j["n_timesteps"] = std::numeric_limits<size_t>::max(); // unlimited
            // else j["n_timesteps"] = j["filenames"].size();
            j["n_timesteps"] = std::numeric_limits<size_t>::max();
          }

          // components
          const size_t nv = j["variables"].size();
          const std::vector<int> components(nv, 1);
          j["components"] = components;
#else
          fatal(FTK_ERR_NOT_BUILT_WITH_ADIOS2);
#endif
        } else if (j["format"] == "bp3") {
#if FTK_HAVE_ADIOS1
          if (missing_variables) 
            fatal("missing variables for bp3");
          const std::string varname0 = j["variables"][0];

          adios_read_init_method( ADIOS_READ_METHOD_BP, comm, "" );
          ADIOS_FILE *fp = adios_read_open_file(filename0.c_str(), ADIOS_READ_METHOD_BP, comm);
  
          ADIOS_VARINFO *avi = adios_inq_var(fp, varname0.c_str());
          adios_inq_var_stat(fp, avi, 0, 0);
          adios_inq_var_blockinfo(fp, avi);
          adios_inq_var_meshinfo(fp, avi);
  
          int nt = 1;
          uint64_t st[4] = {0, 0, 0, 0}, sz[4] = {0, 0, 0, 0};
          std::vector<size_t> mydims;
          
          for (int i = 0; i < avi->ndim; i++) {
            st[i] = 0;
            sz[i] = avi->dims[i];
            nt = nt * sz[i];
            mydims.push_back(sz[i]);
          }
          std::reverse(mydims.begin(), mydims.end());
  
          adios_read_finalize_method (ADIOS_READ_METHOD_BP);
          adios_read_close(fp);
          
          if (j.contains("dimensions"))
            warn("ignoring bp3 dimensions");
          j["dimensions"] = mydims;

          // components
          const size_t nv = j["variables"].size();
          const std::vector<int> components(nv, 1);
          j["components"] = components;
          
          if (j.contains("n_timesteps") && j["n_timesteps"].is_number())
            j["n_timesteps"] = std::min(j["n_timesteps"].template get<size_t>(), j["filenames"].size());
          else 
            j["n_timesteps"] = j["filenames"].size();
#else
          fatal(FTK_ERR_NOT_BUILT_WITH_ADIOS1);
#endif
        } 
      } else fatal("missing filenames");
    } else fatal("invalid input type");
  } else fatal("missing `type'");

  if (j.contains("clamp")) {
    if (j["clamp"].is_array()) {
      auto jc = j["clamp"];
      if (jc.size() == 2) { 
        for (int i = 0; i < 2; i ++) {
          if (jc[0].is_number()) { // ok
          } else fatal("invalid clmap");
        }
        if (jc[0] > jc[1]) 
          fatal("invalid clamp: min is greater than max");
      } else fatal("invalid clamp");
    }
  }

  if (j.contains("perturbation")) {
    if (j["perturbation"].is_number()) { // ok
    } else fatal("invalid perturbation");
  }

  if (j.contains("temporal-smoothing-kernel")) {
    if (j["temporal-smoothing-kernel"].is_number()) {
      if (j.contains("temporal-smoothing-kernel-size")) {
        if (!j["temporal-smoothing-kernel-size"].is_number()) 
          fatal("invalid temporal smoothing kernel size");
      } else 
        j["temporal-smoothing-kernel-size"] = 5; // default value
    } else 
      fatal("invalid temporal smoothing kernel");
  }

  if (j.contains("spatial-smoothing-kernel")) {
    if (j["spatial-smoothing-kernel"].is_number()) {
      if (j.contains("spatial-smoothing-kernel-size")) {
        if (!j["spatial-smoothing-kernel-size"].is_number()) 
          fatal("invalid spatial smoothing kernel size");
      } else 
        j["spatial-smoothing-kernel-size"] = 3; // default value
    } else 
      fatal("invalid spatial smoothing kernel");
  }

  // std::cerr << std::endl 
  //           << j << std::endl;
}

template <typename T>
std::vector<size_t> ndarray_stream<T>::shape() const
{
  std::vector<size_t> shape = j["dimensions"];
  const size_t nc = n_components();
  if (nc > 1)
    shape.insert(shape.begin(), nc);
  // if (is_multi_component()) 
  //   shape.insert(shape.begin(), n_variables());
  return shape;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_file(int k)
{
  const std::string fmt = j["format"];
  if (fmt == "float32") 
    return request_timestep_file_binary<float>(k);
  else if (fmt == "float64")
    return request_timestep_file_binary<double>(k);
  else if (fmt == "vti")
    return request_timestep_file_vti(k);
  else if (fmt == "vtu" || fmt == "vtu_resample" || fmt == "pvtu" || fmt == "pvtu_resample")
    return request_timestep_file_vtu(k);
  else if (fmt == "nc")
    return request_timestep_file_nc(k);
  else if (fmt == "pnc")
    return request_timestep_file_pnc(k);
  else if (fmt == "h5")
    return request_timestep_file_h5(k);
  else if (fmt == "bp3")
    return request_timestep_file_bp3(k);
  else if (fmt == "bp4")
    return request_timestep_file_bp4(k);
  else return ndarray<T>();
}

template <typename T>
template <typename T1>
ndarray<T> ndarray_stream<T>::request_timestep_file_binary(int k)
{
  // TODO FIXME multicomponent binary inputs
  const std::string filename = j["filenames"][k];

  if (comm.size() > 1) { // MPI-IO
    ndarray<T> array;
#if FTK_HAVE_MPI
    MPI_File fp;
    MPI_File_open(comm, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fp);

    std::vector<int> dsz = j["dimensions"], bst, bsz;
    if (!part)
      ext = shape();

    std::vector<size_t> sst = ext.starts(), ssz = ext.sizes();
    for (int i = 0; i < dsz.size(); i ++) {
      bst.push_back(ext.start(i));
      bsz.push_back(ext.size(i));
    }

    std::reverse(dsz.begin(), dsz.end());
    std::reverse(bst.begin(), bst.end());
    std::reverse(bsz.begin(), bsz.end());
    
    ftk::ndarray<T1> array1;
    array1.reshape(ssz);
    // array1.transpose();

    MPI_Datatype dtype;
    MPI_Type_create_subarray(
        dsz.size(), // ndims
        &dsz[0], 
        &bsz[0],
        &bst[0],
        MPI_ORDER_C,
        array1.mpi_datatype(),
        &dtype);
    MPI_Type_commit(&dtype);

    MPI_File_set_view(fp, 0, array1.mpi_datatype(), dtype, "native", MPI_INFO_NULL);

    // MPI_File_read(fp, array1.data(), array1.size(), 
    //     array1.mpi_datatype(), MPI_STATUS_IGNORE);
    MPI_File_read_all(fp, array1.data(), array1.size(), 
        array1.mpi_datatype(), MPI_STATUS_IGNORE);

    MPI_File_close(&fp);
    MPI_Type_free(&dtype);

    array.from_array(array1);
#else
    fatal(FTK_ERR_NOT_BUILT_WITH_MPI);
#endif
    return array;
  } else {
    ftk::ndarray<T> array(shape());

    ftk::ndarray<T1> array1(shape());
    array1.read_binary_file(filename);
    array.from_array(array1);
    
    return array;
  }
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_file_vtu(int k)
{
  ftk::ndarray<T> array;
  const std::string filename = j["filenames"][k];
  bool resample = false;

#if FTK_HAVE_VTK
  if (comm.rank() == 0) {
    vtkSmartPointer<vtkUnstructuredGrid> grid;
    if (j["format" ] == "vtu" || j["format"] == "vtu_resample") {
      if (j["format"] == "vtu_resample")
        resample = true;
      vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
      reader->SetFileName(filename.c_str());
      reader->Update();
      grid = reader->GetOutput();
    } else if (j["format"] == "pvtu" || j["format"] == "pvtu_resample") {
      if (j["format"] == "pvtu_resample")
        resample = true;
      vtkSmartPointer<vtkXMLPUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();
      reader->SetFileName(filename.c_str());
      reader->Update();
      grid = reader->GetOutput();
    }

    if (resample) {
      vtkSmartPointer<vtkResampleToImage> resample = vtkResampleToImage::New();

      std::array<int, 3> dims = {1, 1, 1};
      for (int i = 0; i < j["dimensions"].size(); i ++)
        dims[i] = j["dimensions"][i];
      resample->SetSamplingDimensions(dims[0], dims[1], dims[2]);

      if (j.contains("bounds") && j["bounds"].is_array()) {
        auto jb = j["bounds"];
        std::array<double, 6> b{0};
        for (int i = 0; i < jb.size(); i ++)
          b[i] = jb[i];

        resample->SetUseInputBounds(false);
        resample->SetSamplingBounds(b[0], b[1], b[2], b[3], b[4], b[5]);
      } else {
        resample->SetUseInputBounds(true);
      }
      resample->SetInputDataObject(grid);
      resample->Update();

      vtkSmartPointer<vtkImageData> vti = resample->GetOutput();
      // vti->PrintSelf(std::cerr, vtkIndent(2));

      // read arrays from vti
      const int nv = n_variables();
      if (nv == 1) { //  (is_single_component()) {
        array.from_vtk_image_data(vti, j["variables"][0]);
      } else {
        std::vector<ftk::ndarray<T>> arrays(nv);
        for (int i = 0; i < nv; i ++)
          arrays[i].from_vtk_image_data(vti, j["variables"][i]);

        array.reshape(shape());
        for (int i = 0; i < arrays[0].nelem(); i ++) {
          for (int j = 0; j <nv; j ++) {
            array[i*nv+j] = arrays[j][i];
          }
        }

        array.set_multicomponents();
      }
    } else { // no resample
      const int nv = n_variables();
      if (nv == 1) { //  (is_single_component()) {
        array.from_vtu(grid, j["variables"][0]);
      } else {
        std::vector<ftk::ndarray<T>> arrays(nv);
        for (int i = 0; i < nv; i ++)
          arrays[i].from_vtu(grid, j["variables"][i]);

        array.reshape(shape());
        for (int i = 0; i < arrays[0].nelem(); i ++) {
          for (int j = 0; j <nv; j ++) {
            array[i*nv+j] = arrays[j][i];
          }
        }

        array.set_multicomponents();
      }
    }
  }
  
  if (comm.size() > 1)
    diy::mpi::bcastv(comm, array);
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
  
  if (part)
    return array.subarray(ext);
  else 
    return array;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_file_vti(int k)
{
  ftk::ndarray<T> array;
  const std::string filename = j["filenames"][k];
  // std::cerr << j << std::endl;
#if FTK_HAVE_VTK
  if (comm.rank() == 0) {
    const int nv = n_variables();
    if (nv == 1) { //  (is_single_component()) {
      array.read_vtk_image_data_file(filename, j["variables"][0]);
    } else {
      std::vector<ftk::ndarray<T>> arrays(nv);
      for (int i = 0; i < nv; i ++)
        arrays[i].read_vtk_image_data_file(filename, j["variables"][i]);

      array.reshape(shape());
      for (int i = 0; i < arrays[0].nelem(); i ++) {
        for (int j = 0; j <nv; j ++) {
          array[i*nv+j] = arrays[j][i];
        }
      }

      array.set_multicomponents();
    }
  }
  
  if (comm.size() > 1)
    diy::mpi::bcastv(comm, array);
  
  if (part)
    return array.subarray(ext);
  else 
    return array;

  // std::cerr << array.shape() << ", " << array.multicomponents() << std::endl;
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
  return array;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_file_pnc(int k)
{
  ndarray<T> array;
#if FTK_HAVE_PNETCDF
  fatal("FTK I/O with parallel-netcdf not fully implemented");
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_PNETCDF);
#endif
  return array;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_file_nc(int k)
{
  size_t fid = 0, offset = 0;
  size_t starts[4] = {0}, sizes[4] = {0};

  std::vector<size_t> dsz = j["dimensions"], bst, bsz;

  if (part) {
    bst = ext.starts();
    bsz = ext.sizes();
  } else {
    bst.resize(dsz.size(), 0);
    bsz = dsz;
  }
  
  // determine which file to read
  if (j.contains("first_timestep_per_file")) {
    std::vector<int> first_timestep_per_file = j["first_timestep_per_file"], 
      timesteps_per_file = j["timesteps_per_file"];
    const size_t nf = first_timestep_per_file.size();

    for (size_t i = 0; i < nf; i ++) {
      if (k >= first_timestep_per_file[i] && k < first_timestep_per_file[i] + timesteps_per_file[i]) {
        fid = i;
        offset = k - first_timestep_per_file[i];
        break;
      }
    }

    starts[0] = offset; sizes[0] = 1;
    if (bsz.size() == 1) {
      starts[1] = bst[0];
      sizes[1] = bsz[0];
    } else if (bsz.size() == 2) {
      starts[1] = bst[1];
      starts[2] = bst[0];
      sizes[1] = bsz[1];
      sizes[2] = bsz[0];
    } else { // 3D
      starts[1] = bst[2];
      starts[2] = bst[1];
      starts[3] = bst[0];
      sizes[1] = bsz[2];
      sizes[2] = bsz[1];
      sizes[3] = bsz[0];
    }
  } else {
    fid = k;
    if (bsz.size() == 1) {
      starts[0] = bst[0];
      sizes[0] = bsz[0];
    } else if (bsz.size() == 2) {
      starts[0] = bst[1];
      starts[1] = bst[0];
      sizes[0] = bsz[1];
      sizes[1] = bsz[0];
    } else { // 4D
      starts[0] = bst[2];
      starts[1] = bst[1];
      starts[2] = bst[0];
      sizes[0] = bsz[2];
      sizes[1] = bsz[1];
      sizes[2] = bsz[0];
    }
  }

  // fprintf(stderr, "rank=%d, reading nc, st=%zu, %zu, %zu, %zu, sz=%zu, %zu, %zu, %zu\n", 
  //     comm.rank(),
  //     starts[0], starts[1], starts[2], starts[3], 
  //     sizes[0], sizes[1], sizes[2], sizes[3]);
  // fprintf(stderr, "k=%d, offset=%zu, fid=%zu\n", k, offset, fid);

  ftk::ndarray<T> array;
  const std::string filename = j["filenames"][fid];
#if FTK_HAVE_NETCDF
  const int nv = n_variables();
  if (nv == 1) { // all data in one single variable; channels are automatically handled in ndarray
    array.read_netcdf(filename, j["variables"][0], starts, sizes);
    // array.reshape(shape()); // ncdims may not be equal to nd
    array.reshape(bsz); // ncdims may not be equal to nd
  } else { // u, v, w in separate variables
    const int nv = n_variables();
    std::vector<ftk::ndarray<T>> arrays(nv);
    for (int i = 0; i < nv; i ++) {
      arrays[i].read_netcdf(filename, j["variables"][i], starts, sizes);
      // arrays[i].transpose(); // TODO FIXME: this is strange.. works for MPAS-O; will see if works for other applications
    }

    // array.reshape(shape());
    bsz.insert(bsz.begin(), nv);
    array.reshape(bsz); 
    for (int i = 0; i < arrays[0].nelem(); i ++) {
      for (int j = 0; j <nv; j ++) {
        array[i*nv+j] = arrays[j][i];
      }
    }

    array.set_multicomponents();
  }
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
  return array;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_file_h5(int k)
{
  ftk::ndarray<T> array;
  const std::string filename = j["filenames"][k];
#if FTK_HAVE_HDF5
  const int nc = n_components();
  if (nc == 1) { // all data in one single-component variable; channels are automatically handled in ndarray
    array.read_h5(filename, j["variables"][0]);
    array.reshape(shape()); // ncdims may not be equal to nd
  } else { // u, v, w in separate variables
    const int nv = n_variables();
    std::vector<ftk::ndarray<T>> arrays(nv);
    for (int i = 0; i < nv; i ++)
      arrays[i].read_h5(filename, j["variables"][i]);

    array = ndarray<T>::concat(arrays);
    array.set_multicomponents();
  }
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_HDF5);
#endif
  return array;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_file_bp3(int k)
{
  ftk::ndarray<T> array;
  const std::string filename = j["filenames"][k];
#if FTK_HAVE_ADIOS1
  const int nc = n_components();
  if (nc == 1) { // all data in one single-component variable; channels are automatically handled in ndarray
    array.read_bp_legacy(filename, j["variables"][0], comm);
  } else { // u, v, w in separate variables
    const int nv = n_variables();
    std::vector<ftk::ndarray<T>> arrays(nv);
    for (int i = 0; i < nv; i ++)
      arrays[i].read_bp_legacy(filename, j["variables"][i], comm);

    array = ndarray<T>::concat(arrays);
    array.set_multicomponents();
  }
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_ADIOS1);
#endif
  return array;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_file_bp4(int k)
{
  ftk::ndarray<T> array;
#if FTK_HAVE_ADIOS2
  fprintf(stderr, "requesting bp4 timestep %d\n", k);
  while (k != 0) { // the first step was already beginned with configuration
    auto status = adios_reader.BeginStep( adios2::StepMode::Read, 10.f );
    // std::cerr << "status=" << status << std::endl;
    // fprintf(stderr, "status=%d\n", status);
    if (status == adios2::StepStatus::NotReady) {
      usleep(100000);
      continue;
    } else if (status == adios2::StepStatus::OK) {
      break;
    } else if (status == adios2::StepStatus::EndOfStream) {
      fprintf(stderr, "end of stream, current_file_id=%d\n", current_file_id);
      current_file_id ++; 
      if (current_file_id < j["filenames"].size()) {
        const std::string filename = j["filenames"][current_file_id];
        adios_reader.Close();
        adios_reader = adios_io.Open(filename, adios2::Mode::Read);
        continue;
      } else {
        return array; // EOF
        break;
      }
    } else { // other error
      fatal("adios2 error");
    }
  }

  const int nc = n_components();
  if (nc == 1) { // all data in one single-component variable; channels are automatically handled in ndarray
    // array.read_bp(filename, j["variables"][0], comm);
    auto var = adios_io.template InquireVariable<T>(j["variables"][0]);
    array.read_bp(adios_io, adios_reader, var); // adios_vars[0]);
  } else { // u, v, w in separate variables
    const int nv = n_variables();
    std::vector<ftk::ndarray<T>> arrays(nv);
    for (int i = 0; i < nv; i ++) {
      // arrays[i].read_bp(filename, j["variables"][i], comm);
      auto var = adios_io.template InquireVariable<T>(j["variables"][0]);
      array.read_bp(adios_io, adios_reader, var); // adios_vars[i]);
    }

    array = ndarray<T>::concat(arrays);
    array.set_multicomponents();
  }

  adios_reader.EndStep();
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_ADIOS2);
#endif
  return array;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic(int k)
{
  if (j["name"] == "woven") 
    return request_timestep_synthetic_woven(k);
  else if (j["name"] == "moving_extremum_2d")
    return request_timestep_synthetic_moving_extremum_2d(k);
  else if (j["name"] == "moving_extremum_3d")
    return request_timestep_synthetic_moving_extremum_3d(k);
  else if (j["name"] == "moving_ramp_3d")
    return request_timestep_synthetic_moving_ramp_3d(k);
  else if (j["name"] == "moving_dual_ramp_3d")
    return request_timestep_synthetic_moving_dual_ramp_3d(k);
  else if (j["name"] == "double_gyre")
    return request_timestep_synthetic_double_gyre(k);
  else if (j["name"] == "merger_2d")
    return request_timestep_synthetic_merger_2d(k);
  else if (j["name"] == "volcano_2d")
    return request_timestep_synthetic_volcano_2d(k);
  else if (j["name"] == "tornado")
    return request_timestep_synthetic_tornado(k);
  return ndarray<T>();
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_woven(int k) 
{
  const int nt = j["n_timesteps"];
  const T t = nt == 1 ? 0.0 : double(k)/(nt-1);
  const size_t DW = j["dimensions"][0], DH = j["dimensions"][1];

  if (part) {
    lattice domain({0, 0}, {DW, DH});
    return ftk::synthetic_woven_2D_part<T>(domain, ext, t);
  } else {
    return ftk::synthetic_woven_2D<T>(DW, DH, t);
  }
  // if (part) return arr.subarray(ext);
  // else return arr;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_moving_extremum_2d(int k) 
{
  const std::vector<size_t> shape({j["dimensions"][0], j["dimensions"][1]});
  const T x0[2] = {j["x0"][0], j["x0"][1]}, 
          dir[2] = {j["dir"][0], j["dir"][1]};

  return ftk::synthetic_moving_extremum<T, 2>(shape, x0, dir, T(k));
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_volcano_2d(int k) 
{
  const std::vector<size_t> shape({j["dimensions"][0], j["dimensions"][1]});
  const T x0[2] = {j["x0"][0], j["x0"][1]}, 
          radius = j["radius"];

  return ftk::synthetic_volcano<T, 2>(shape, x0, radius+k);
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_moving_extremum_3d(int k) 
{
  const std::vector<size_t> shape({j["dimensions"][0], j["dimensions"][1], j["dimensions"][2]});
  const T x0[3] = {j["x0"][0], j["x0"][1], j["x0"][2]}, 
          dir[3] = {j["dir"][0], j["dir"][1], j["dir"][2]};

  return ftk::synthetic_moving_extremum<T, 3>(shape, x0, dir, T(k));
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_moving_ramp_3d(int k) 
{
  const std::vector<size_t> shape({j["dimensions"][0], j["dimensions"][1], j["dimensions"][2]});
  const T x0 = j["x0"], rate = j["rate"];

  return ftk::synthetic_moving_ramp<T, 3>(shape, x0, rate, T(k));
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_moving_dual_ramp_3d(int k) 
{
  const std::vector<size_t> shape({j["dimensions"][0], j["dimensions"][1], j["dimensions"][2]});
  const T x0 = j["x0"], rate = j["rate"], offset = j["offset"];

  return ftk::synthetic_moving_dual_ramp<T, 3>(shape, x0, rate, T(k) + offset);
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_merger_2d(int k) 
{
  return synthetic_merger_2D(j["dimensions"][0], j["dimensions"][1], T(k)*0.1);
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_double_gyre(int k) 
{
  // TODO: allowing parameter configuration
  const double time_scale = j["time_scale"];
  const T A = 0.1, Omega = M_PI * 2, Eps = 0.25;
  const T time = k * time_scale; 
  const size_t DW = j["dimensions"][0], 
               DH = j["dimensions"][1];
 
  auto result = synthetic_double_gyre<T>(DW, DH, time, A, Omega, Eps);
  // std::cerr << result.multicomponents() << std::endl;
  return result;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_tornado(int k) 
{
  return ftk::synthetic_tornado<T>(
      j["dimensions"][0],
      j["dimensions"][1],
      j["dimensions"][2], 
      k);
}

template <typename T>
void ndarray_stream<T>::modified_callback(int k, const ndarray<T> &array)
{
  auto f = [&](const ndarray<T> &array) {
    if (j.contains("temporal-smoothing-kernel"))
      temporal_filter.push(array);
    else
      callback(k, array);
  };
   
  auto f1 = [&](const ndarray<T> &array) { // handling clamp
    if (j.contains("clamp")) {
      auto array1 = array;
      array1.clamp(j["clamp"][0], j["clamp"][1]);
      f(array1);
    } else 
      f(array);
  };

  auto f2 = [&](const ndarray<T> &array) { // handling perturbation
    if (j.contains("perturbation")) { 
      auto array1 = array; 
      array1.perturb(j["perturbation"]);
      f1(array1);
    } else 
      f1(array);
  };

  if (j.contains("spatial-smoothing-kernel")) {
    const int ksize = j["spatial-smoothing-kernel-size"];
    const T sigma = j["spatial-smoothing-kernel"];
    ndarray<T> array1 = conv_gaussian(array, sigma, ksize, ksize/2);
    f2(array1);
  } else
    f2(array);
}

template <typename T>
void ndarray_stream<T>::start()
{
  if (!callback) 
    fatal("callback function not set");

  if (j.contains("temporal-smoothing-kernel")) {
    temporal_filter.set_gaussian_kernel(j["temporal-smoothing-kernel"], j["temporal-smoothing-kernel-size"]);
    temporal_filter.set_callback(callback);
  }

  const int nt = j["n_timesteps"];
  const bool async = j.contains("async");

  if (async) {
    std::condition_variable cond;
    std::mutex mtx;
    bool all_done = false;
    int current_timestep;
    ndarray<T> array;

    std::thread producer([&]() {
      for (int i = 0; i < nt; i ++) {
        // fprintf(stderr, "reading timestep %d\n", i);
        ndarray<T> array1;

        auto t1 = std::chrono::high_resolution_clock::now();
        if (j["type"] == "synthetic") {
          array1 = request_timestep_synthetic(i);
        } else if (j["type"] == "file") {
          array1 = request_timestep_file(i);
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        float t = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count() * 1e-9;
        fprintf(stderr, "timestep=%d, t_io=%f\n", i, t);

        {
          std::unique_lock<std::mutex> lock(mtx);

          array = array1;
          current_timestep = i;
          if (i == nt-1 || array1.empty()) {
            all_done = true;
          }
          cond.notify_one();
        }
      }
      // fprintf(stderr, "exiting reader thread..\n");
    });

    { // consumer
      std::unique_lock<std::mutex> lock(mtx);
      while (!all_done) {
        cond.wait(lock); // , [&](){return all_done;});
      
        auto t1 = std::chrono::high_resolution_clock::now();
        modified_callback(current_timestep, array);
        auto t2 = std::chrono::high_resolution_clock::now();
        float t = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count() * 1e-9;
        fprintf(stderr, "timestep=%d, t_compute=%f\n", current_timestep, t);
      }
    }

    producer.join();
  } else {
    for (size_t i = 0; i < nt; i ++) {
      auto t0 = std::chrono::high_resolution_clock::now();
      ndarray<T> array;
      if (j["type"] == "synthetic") {
        array = request_timestep_synthetic(i);
        // fprintf(stderr, "requested data ncd=%zu, ncd1=%zu\n", array.multicomponents(), array1.multicomponents());
      }
      else if (j["type"] == "file") {
        array = request_timestep_file(i);
        if (array.empty()) {
          fprintf(stderr, "got empty array; all files are read.\n"); 
          break;
        }
      }
      auto t1 = std::chrono::high_resolution_clock::now();

      modified_callback(i, array);
      auto t2 = std::chrono::high_resolution_clock::now();

      float t_io = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() * 1e-9,
            t_compute = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() * 1e-9;
      
      fprintf(stderr, "timestep=%zu, t_io=%f, t_compute=%f\n", i, t_io, t_compute);
    }
  }
}
 
template <typename T>
void ndarray_stream<T>::finish() 
{
  if (j.contains("temporal-smoothing-kernel"))
    temporal_filter.finish();
}

template <typename T>
bool ndarray_stream<T>::is_partial_read_supported() const
{
  if (j["type"] == "synthetic") {
    if (j["name"] == "woven") return true;
    else return false;
  } else { // file
    const std::string fmt = j["format"];
    if (fmt == "float32" || fmt == "float64") 
      return true; 
    else if (fmt == "vti")
      return true; // compromised
    else if (fmt == "vtu")
      return true; // compromised
    else if (fmt == "nc") 
      return true; // performance depending on if parallel IO is supported
    else if (fmt == "pnc")
      return false; // TODO
    else if (fmt == "h5")
      return false; // TODO
    else if (fmt == "bp3")
      return false; // TODO
    else if (fmt == "bp4")
      return false; // TODO
    else 
      return false;
  } 
}

} // namespace ftk

#endif
