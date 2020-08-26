#ifndef _FTK_REGULAR_ARRAY_STREAM_HH
#define _FTK_REGULAR_ARRAY_STREAM_HH

#include <ftk/object.hh>
#include <fstream>
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
  //  - format (required if type is file and format is float32/float64), string.  If not 
  //    given, the format will be determined by the filename extension.  The value of this 
  //    field must be one of the follows: vti, nc, h5, float32, float64.
  //  - variables (required if format is nc/h5, optional for vti), array of strings.
  //    - the number of components is the length of the array.
  //    - if not given, the defaulat value is ["scalar"]
  //  - dimensions (required if format is floaot32/float64), array of integers, e.g.  
  //    [width, height, depth]
  //  - n_timesteps, integer.  The default is 32 for synthetic data; the number can be 
  //    automatically derived from file; the number can be automatically derived from files

  void set_input_source_json_file(const std::string& filename);
  void set_input_source_json(const json& j_);
  const json& get_json() const {return j;}
  
  void start();
  void finish();

  void set_callback(std::function<void(int, const ndarray<T>&)> f) {callback = f;}

  bool is_single_component() const { return n_components() == 1; }
  bool is_multi_component() const { return n_components() > 1; }
  size_t n_components() const { return j["variables"].size(); }
  size_t n_dimensions() const {
    // if (j.contains("nc_has_unlimited_time_dimension")) return j["dimensions"].size() - 1;
    // else return j["dimensions"].size(); 
    return j["dimensions"].size(); 
  }

  std::vector<size_t> shape() const;

protected:
  ndarray<T> request_timestep_file(int k);
  ndarray<T> request_timestep_file_nc(int k);
  ndarray<T> request_timestep_file_vti(int k);
  ndarray<T> request_timestep_file_h5(int k);
  template <typename T1> ndarray<T> request_timestep_file_binary(int k);

  ndarray<T> request_timestep_synthetic(int k);
  ndarray<T> request_timestep_synthetic_woven(int k);
  ndarray<T> request_timestep_synthetic_moving_extremum_2d(int k);
  ndarray<T> request_timestep_synthetic_moving_extremum_3d(int k);
  ndarray<T> request_timestep_synthetic_double_gyre(int k);
  ndarray<T> request_timestep_synthetic_merger_2d(int k);
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
            j["x0"] = {0.5, 0.5, 0.5};

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
            j["dir"] = {0.1, 0.1, 0.1};
          
        } else if (j["name"] == "double_gyre") {
          default_nd = 2;
        } else if (j["name"] == "merger_2d") {
          j["variables"] = {"scalar"};
          default_nd = 2;
          default_dims[0] = 32; 
          default_dims[1] = 32;
          default_n_timesteps = 100;
        } else if (j["name"] == "tornado") {
          default_nd = 3;
          j["variables"] = {"u", "v", "w"};
        } else fatal("synthetic case not available.");
      } else fatal("synthetic case name not given.");
     
      if (missing_dimensions) {
        std::vector<int> dims;
        for (int i = 0; i < default_nd; i ++)
          dims.push_back(default_dims[i]);
        j["dimensions"] = dims;
      }

      // if (j.contains("n_timesteps")) assert(j["n_timesteps"] != 0);
      if (!j.contains("n_timesteps"))
        j["n_timesteps"] = default_n_timesteps;
    } else if (j["type"] == "file") {
      if (j.contains("filenames")) {
        if (!j["filenames"].is_array()) {
          auto filenames = ftk::ndarray<double>::glob(j["filenames"]);
          if (filenames.empty()) fatal("unable to find matching filename(s).");
          if (j.contains("n_timesteps")) filenames.resize(j["n_timesteps"]);
          else j["n_timesteps"] = filenames.size();
          j["filenames"] = filenames;
        }
        const std::string filename0 = j["filenames"][0];

        if (!j.contains("format")) { // probing file format
          if (ends_with(filename0, "vti")) j["format"] = "vti";
          else if (ends_with(filename0, "nc")) j["format"] = "nc";
          else if (ends_with(filename0, "h5")) j["format"] = "h5";
          else fatal("unable to determine file format.");
        }

        if (j["format"] == "float32" || j["format"] == "float64") {
          if (missing_dimensions) 
            fatal("missing dimensions.");
        
          if (missing_variables) { // add one single component in the array
            const std::string var("scalar");
            const std::vector<std::string> vars = {var};
            j["variables"] = vars;
          }
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
          
          // determine number timesteps
          if (j.contains("n_timesteps") && j["n_timesteps"].is_number())
            j["n_timesteps"] = std::min(j["n_timesteps"].template get<size_t>(), j["filenames"].size());
          else 
            j["n_timesteps"] = j["filenames"].size();
#else
          fatal("FTK not compiled with VTK.");
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

          // determin spatial dimensions
          int dimids[4]; 
          size_t dimlens[4];
          NC_SAFE_CALL( nc_inq_vardimid(ncid, varids[0], dimids) );
          for (int i = 0; i < ncdims; i ++) 
            NC_SAFE_CALL( nc_inq_dimlen(ncid, dimids[i], &dimlens[i]) );

          int unlimited_recid;
          NC_SAFE_CALL( nc_inq_unlimdim(ncid, &unlimited_recid) );
          NC_SAFE_CALL( nc_close(ncid) );
          if (unlimited_recid >= 0)
            j["nc_has_unlimited_time_dimension"] = true;
          
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

          // determine number timesteps
          if (j.contains("n_timesteps") && j["n_timesteps"].is_number())
            j["n_timesteps"] = std::min(j["n_timesteps"].template get<size_t>(), j["filenames"].size());
          else 
            j["n_timesteps"] = j["filenames"].size();
#else
          fatal("FTK not compiled with NetCDF.");
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
#else
          fatal("FTK not compiled with HDF5.");
          // fatal("array stream w/ h5 not implemented yet");
#endif
        }
      } else fatal("missing filenames");
    } else fatal("invalid input type");
  } else fatal("missing `type'");
 

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
  if (is_multi_component()) 
    shape.insert(shape.begin(), n_components());
  return shape;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_file(int k)
{
  if (j["format"] == "float32") 
    return request_timestep_file_binary<float>(k);
  else if (j["format"] == "float64")
    return request_timestep_file_binary<double>(k);
  else if (j["format"] == "vti")
    return request_timestep_file_vti(k);
  else if (j["format"] == "nc")
    return request_timestep_file_nc(k);
  else if (j["format"] == "h5")
    return request_timestep_file_h5(k);
  else return ndarray<T>();
}

template <typename T>
template <typename T1>
ndarray<T> ndarray_stream<T>::request_timestep_file_binary(int k)
{
  const std::string filename = j["filenames"][k];
  ftk::ndarray<T1> array1(shape());
  array1.from_binary_file(filename);
  
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
  const int nv = n_components();
  if (is_single_component()) {
    array.from_vtk_image_data_file(filename, j["variables"][0]);
  } else {
    std::vector<ftk::ndarray<T>> arrays(nv);
    for (int i = 0; i < nv; i ++)
      arrays[i].from_vtk_image_data_file(filename, j["variables"][i]);

    array.reshape(shape());
    for (int i = 0; i < arrays[0].nelem(); i ++) {
      for (int j = 0; j <nv; j ++) {
        array[i*nv+j] = arrays[j][i];
      }
    }
  }
#else
  fatal("FTK not compiled with VTK.");
#endif
  return array;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_file_nc(int k)
{
  ftk::ndarray<T> array;
  const std::string filename = j["filenames"][k];
#if FTK_HAVE_NETCDF
  if (is_single_component()) { // all data in one single variable; channels are automatically handled in ndarray
    array.from_netcdf(filename, j["variables"][0]);
    array.reshape(shape()); // ncdims may not be equal to nd
  } else { // u, v, w in separate variables
    const int nv = n_components();
    std::vector<ftk::ndarray<T>> arrays(nv);
    for (int i = 0; i < nv; i ++)
      arrays[i].from_netcdf(filename, j["variables"][i]);

    array.reshape(shape());
    for (int i = 0; i < arrays[0].nelem(); i ++) {
      for (int j = 0; j <nv; j ++) {
        array[i*nv+j] = arrays[j][i];
      }
    }
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
  if (is_single_component()) { // all data in one single variable; channels are automatically handled in ndarray
    array.from_h5(filename, j["variables"][0]);
    array.reshape(shape()); // ncdims may not be equal to nd
  } else { // u, v, w in separate variables
    const int nv = n_components();
    std::vector<ftk::ndarray<T>> arrays(nv);
    for (int i = 0; i < nv; i ++)
      arrays[i].from_h5(filename, j["variables"][i]);

    array = ndarray<T>::concat(arrays);
  }
#else
  fatal("FTK not compiled with HDF5");
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
  else if (j["name"] == "double_gyre")
    return request_timestep_synthetic_double_gyre(k);
  else if (j["name"] == "merger_2d")
    return request_timestep_synthetic_merger_2d(k);
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
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_moving_extremum_3d(int k) 
{
  const std::vector<size_t> shape({j["dimensions"][0], j["dimensions"][1], j["dimensions"][2]});
  const T x0[3] = {j["x0"][0], j["x0"][1], j["x0"][2]}, 
          dir[3] = {j["dir"][0], j["dir"][1], j["dir"][2]};

  return ftk::synthetic_moving_extremum<T, 3>(shape, x0, dir, T(k));
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_merger_2d(int k) 
{
  return synthetic_merger_2D(j["dimensions"][0], j["dimensions"][1], T(k)*0.1);
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_double_gyre(int k) 
{
  return ndarray<T>(); // TODO
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

  if (j.contains("spatial-smoothing-kernel")) {
    const int ksize = j["spatial-smoothing-kernel-size"];
    const T sigma = j["spatial-smoothing-kernel"];
    ndarray<T> array1 = conv_gaussian(array, sigma, ksize, ksize/2);
    f(array1);
  } else
    f(array);
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
    ndarray<T> array;
    if (j["type"] == "synthetic") 
      array = request_timestep_synthetic(i);
    else if (j["type"] == "file")
      array = request_timestep_file(i);
    modified_callback(i, array);
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
