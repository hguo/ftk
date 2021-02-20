#ifndef _FTK_REGULAR_ARRAY_STREAM_HH
#define _FTK_REGULAR_ARRAY_STREAM_HH

#include <ftk/object.hh>
#include <fstream>
#include <chrono>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/filters/streaming_filter.hh>
#include <ftk/external/json.hh>

namespace ftk {
using nlohmann::json;

template <typename T=double>
struct ndarray_stream : public object {
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
  //    field must be one of the follows: vti, nc, h5, float32, float64.
  //  - variables (required if format is nc/h5, optional for vti), array of strings.
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
  ndarray<T> request_timestep_file_vti(int k);
  ndarray<T> request_timestep_file_h5(int k);
  ndarray<T> request_timestep_file_bp(int k);
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

protected:
  json j; // configs, metadata, and everything

  std::function<void(int, const ndarray<T>&)> callback;

  streaming_filter<ndarray<T>, T> temporal_filter;

private:
  static bool ends_with(std::string const & value, std::string const & ending)
  {
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
  }
};

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
        } else {
          auto filenames = glob(j["filenames"]);
          if (filenames.empty()) fatal("unable to find matching filename(s).");
          // if (j.contains("n_timesteps")) filenames.resize(j["n_timesteps"]);
          // else j["n_timesteps"] = filenames.size();
          j["filenames"] = filenames;
        }
        const std::string filename0 = j["filenames"][0];

        if (!j.contains("format")) { // probing file format
          const auto ext = file_extension(filename0);
          
          if (ext == FILE_EXT_VTI) j["format"] = "vti";
          else if (ext == FILE_EXT_NETCDF) j["format"] = "nc";
          else if (ext == FILE_EXT_HDF5) j["format"] = "h5";
          else if (ext == FILE_EXT_BP) j["format"] = "bp";
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
          
          j["n_timesteps"] = j["filenames"].size(); // TODO: we are assuming #timesteps = #filenames
        } else if (j["format"] == "vti") {
#if FTK_HAVE_VTK
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
#else
          fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
        } else if (j["format"] == "nc") {
#if FTK_HAVE_NETCDF
          if (missing_variables) 
            fatal("missing nc variable");

          const int nv = j["variables"].size();
          int ncid, ncdims, nd, varids[nv];

          NC_SAFE_CALL( nc_open(filename0.c_str(), NC_NOWRITE, &ncid) );
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
          
          if (j.contains("dimensions"))
            warn("ignorning dimensions");

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

#else
          fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
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
          std::reverse(dims.begin(), dims.end()); // in h5, the last listed dimension is the fastest-changing dimension
          
          if (j.contains("dimensions"))
            warn("ignoring dimensions");
          j["dimensions"] = dims;

          H5Fclose(fid);
          
          j["n_timesteps"] = j["filenames"].size(); // TODO: we are assuming #timesteps = #filenames
          
          // components
          const size_t nv = j["variables"].size();
          const std::vector<int> components(nv, 1);
          j["components"] = components;
#else
          fatal(FTK_ERR_NOT_BUILT_WITH_HDF5);
#endif
        } else if (j["format"] == "bp") {
#if FTK_HAVE_ADIOS2
          if (missing_variables)
            fatal("missing variables for bp");
          const std::string varname0 = j["variables"][0];

          adios2::ADIOS adios(comm);
          adios2::IO io = adios.DeclareIO("BPReader");
          adios2::Engine reader = io.Open(filename0, adios2::Mode::Read);

          auto var = io.template InquireVariable<T>(varname0);
          if (!var) fatal("variable not found in bp");

          std::vector<size_t> dims(var.Shape());
          if (j.contains("dimensions"))
            warn("ignoring bp dimensions");
          j["dimensions"] = dims;

          reader.Close();

          j["n_timesteps"] = j["filenames"].size();

          // components
          const size_t nv = j["variables"].size();
          const std::vector<int> components(nv, 1);
          j["components"] = components;
#else
          fatal(FTK_ERR_NOT_BUILT_WITH_ADIOS2);
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
  else if (fmt == "nc")
    return request_timestep_file_nc(k);
  else if (fmt == "h5")
    return request_timestep_file_h5(k);
  else if (fmt == "bp")
    return request_timestep_file_bp(k);
  else return ndarray<T>();
}

template <typename T>
template <typename T1>
ndarray<T> ndarray_stream<T>::request_timestep_file_binary(int k)
{
  const std::string filename = j["filenames"][k];
  ftk::ndarray<T1> array1(shape());
  array1.read_binary_file(filename);
  
  ftk::ndarray<T> array(shape());
  array.from_array(array1);

  return array;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_file_vti(int k)
{
  ftk::ndarray<T> array;
  const std::string filename = j["filenames"][k];
  // std::cerr << j << std::endl;
#if FTK_HAVE_VTK
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

  // std::cerr << array.shape() << ", " << array.multicomponents() << std::endl;
#else
  fatal("FTK not compiled with VTK.");
#endif
  return array;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_file_nc(int k)
{
  size_t fid = 0, offset = 0;
  size_t starts[4] = {0}, sizes[4] = {0};
  std::vector<size_t> dims = j["dimensions"];
  
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
    if (dims.size() == 2) {
      sizes[1] = dims[1];
      sizes[2] = dims[0];
    } else {
      sizes[1] = dims[2];
      sizes[2] = dims[1];
      sizes[3] = dims[0];
    }
  } else {
    fid = k;
    if (dims.size() == 2) {
      sizes[0] = dims[1];
      sizes[1] = dims[0];
    } else {
      sizes[0] = dims[2];
      sizes[1] = dims[1];
      sizes[2] = dims[0];
    }
  }

  fprintf(stderr, "st=%zu, %zu, %zu, %zu, sz=%zu, %zu, %zu, %zu\n", 
      starts[0], starts[1], starts[2], starts[3], 
      sizes[0], sizes[1], sizes[2], sizes[3]);
  fprintf(stderr, "k=%d, offset=%zu, fid=%zu\n", k, offset, fid);

  ftk::ndarray<T> array;
  const std::string filename = j["filenames"][fid];
#if FTK_HAVE_NETCDF
  const int nv = n_variables();
  if (nv == 1) { // all data in one single variable; channels are automatically handled in ndarray
    array.read_netcdf(filename, j["variables"][0], starts, sizes);
    array.reshape(shape()); // ncdims may not be equal to nd
  } else { // u, v, w in separate variables
    const int nv = n_variables();
    std::vector<ftk::ndarray<T>> arrays(nv);
    for (int i = 0; i < nv; i ++)
      arrays[i].read_netcdf(filename, j["variables"][i], starts, sizes);

    array.reshape(shape());
    for (int i = 0; i < arrays[0].nelem(); i ++) {
      for (int j = 0; j <nv; j ++) {
        array[i*nv+j] = arrays[j][i];
      }
    }

    array.set_multicomponents();
  }
#else
  fatal("FTK not compiled with netcdf");
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
ndarray<T> ndarray_stream<T>::request_timestep_file_bp(int k)
{
  ftk::ndarray<T> array;
  const std::string filename = j["filenames"][k];
#if FTK_HAVE_ADIOS2
  const int nc = n_components();
  if (nc == 1) { // all data in one single-component variable; channels are automatically handled in ndarray
    array.read_bp(filename, j["variables"][0], comm);
  } else { // u, v, w in separate variables
    const int nv = n_variables();
    std::vector<ftk::ndarray<T>> arrays(nv);
    for (int i = 0; i < nv; i ++)
      arrays[i].read_bp(filename, j["variables"][i]);

    array = ndarray<T>::concat(arrays);
    array.set_multicomponents();
  }
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
  return ftk::synthetic_woven_2D<T>(j["dimensions"][0], j["dimensions"][1], t);
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

  for (int i = 0; i < j["n_timesteps"]; i ++) {
    auto t0 = std::chrono::high_resolution_clock::now();
    ndarray<T> array;
    if (j["type"] == "synthetic") {
      array = request_timestep_synthetic(i);
      // fprintf(stderr, "requested data ncd=%zu, ncd1=%zu\n", array.multicomponents(), array1.multicomponents());
    }
    else if (j["type"] == "file")
      array = request_timestep_file(i);
    auto t1 = std::chrono::high_resolution_clock::now();

    modified_callback(i, array);
    auto t2 = std::chrono::high_resolution_clock::now();

    float t_io = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() * 1e-9,
          t_compute = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() * 1e-9;
    
    fprintf(stderr, "timestep=%d, t_io=%f, t_compute=%f\n", i, t_io, t_compute);
  }
}
 
template <typename T>
void ndarray_stream<T>::finish() 
{
  if (j.contains("temporal-smoothing-kernel"))
    temporal_filter.finish();
}

}

#endif
