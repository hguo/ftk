#ifndef _FTK_REGULAR_ARRAY_STREAM_HH
#define _FTK_REGULAR_ARRAY_STREAM_HH

#include <ftk/ndarray.hh>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/filters/streaming_filter.hh>
#include <ftk/external/json.hh>

namespace ftk {
using nlohmann::json;

template <typename T=double>
struct ndarray_stream {
  void set_input_source_json_file(const std::string& filename);
  void set_input_source_json(const json& j_);
  const json& get_json() const {return j;}
  
  void start();
  void finish();

  void set_callback(std::function<void(int, ndarray<T>&)> f) {callback = f;}

  bool is_single_component() const { return j["variable"].is_string(); }
  bool is_multi_component() const { return j["variable"].is_array(); }
  size_t n_components() const { 
    if (j["variable"].is_string()) return 1;
    else if (j["variable"].is_array()) return j["variable"].size();
    else { fatal("invalid variable"); return 0;}
  }

  std::vector<size_t> shape() const;

protected:
  ndarray<T> request_timestep_file(int k);
  ndarray<T> request_timestep_file_nc(int k);
  ndarray<T> request_timestep_file_vti(int k);
  ndarray<T> request_timestep_file_h5(int k);
  template <typename T1> ndarray<T> request_timestep_binary(int k);

  ndarray<T> request_timestep_synthetic(int k);
  ndarray<T> request_timestep_synthetic_woven(int k);
  ndarray<T> request_timestep_synthetic_double_gyre(int k);
  ndarray<T> request_timestep_synthetic_merger(int k);

protected:
  static void fatal(const std::string& str) {
    std::cerr << "FATAL: " << str << std::endl;
    exit(1);
  }

  static void warn(const std::string& str) {
    std::cerr << "WARN: " << str << std::endl;
  }

protected:
  json j; // configs, metadata, and everything

  int current_timestep = 0;
  std::function<void(int, ndarray<T>&)> callback;

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
}

template <typename T>
void ndarray_stream<T>::set_input_source_json(const json& j_)
{
  j = j_;
  // std::cerr << j << std::endl;
  
  if (j.contains("type")) {
    if (j["type"] == "synthetic") {
      if (j.contains("name")) {
        if (j["name"] == "woven") {
          j["nd"] = 2;
          if (!j.contains("scaling_factor")) j["scalaring_factor"] = 15.0;
        } else if (j["name"] == "double_gyre") {
          j["nd"] = 2;
        } else if (j["name"] == "merger") {
          j["nd"] = 2;
        } else fatal("synthetic case not available.");
      } else fatal("synthetic case name not given.");
     
      // default dimensions
      if (j.contains("width")) assert(j["width"] != 0);
      else j["width"] = 32;
      
      if (j.contains("height")) assert(j["height"] != 0);
      else j["height"] = 32;
        
      if (j["nd"] == 3) {
        if (j.contains("depth")) assert(j["depth"] != 0);
        else j["depth"] = 32;
      }

      // if (j.contains("n_timesteps")) assert(j["n_timesteps"] != 0);
      if (!j.contains("n_timesteps"))
        j["n_timesteps"] = 32;
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

        if (j.contains("variable")) { 
          if (j["variable"].is_array()) { // multicomponent
            for (const auto &v : j["variable"]) {
              if (!v.is_string()) fatal("invalid variable name");
            }
          } else if (j["variable"].is_string()) { // single-component
            
          } else fatal("invalid variable");
        } // else fatal("missing variable");

        if (j["format"] == "float32" || j["format"] == "float64") {
          if (j.contains("nd")) {
            if (j["nd"] != 2 && j["nd"] != 3) fatal("unsupported spatial dimensionality");
          } else fatal("unable to determine spatial dimensionality");

          if ((j["nd"] == 2 && ((!j.contains("width") || !j.contains("height")))) || 
              (j["nd"] == 3 && ((!j.contains("width") || !j.contains("height") || !j.contains("depth")))))
            fatal("width, height, and/or depth not specified.");
        } else if (j["format"] == "vti") {
#if FTK_HAVE_VTK
          vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
          reader->SetFileName(filename0.c_str());
          reader->Update();

          vtkSmartPointer<vtkImageData> image = reader->GetOutput();

          // determine dimensionality
          int imageNd = image->GetDataDimension();
          if (j.contains("nd") && j["nd"].is_number()) {
            if (j["nd"] != imageNd)
              fatal("Data dimensionality and given dimensionality mismatch.");
          } else {
            j["nd"] = imageNd;
            if (!(j["nd"] == 2 || j["nd"] == 3))
              fatal("Unsupported VTI data dimensionality.");
          }

          if (j.contains("width") || j.contains("height") || j.contains("depth")) 
            warn("Given data dimensions are ignored.");
          j["width"] = image->GetDimensions()[0];
          j["height"] = image->GetDimensions()[1];
          if (j["nd"] == 3)
            j["depth"] = image->GetDimensions()[2];

          if (!j.contains("variable")) {
            j["variable"] = image->GetPointData()->GetArrayName(0);
          } else if (j["variable"].is_string()) {
            int varid;
            image->GetPointData()->GetArray(j["variable"].template get<std::string>().c_str(), varid);
          } else if (j["variable"].is_array()) {
            for (int i = 0; i < j["variable"].size(); i ++) {
              int varid;
              image->GetPointData()->GetArray(j["variable"][i].template get<std::string>().c_str(), varid);
              if (varid < 0) fatal("cannot find variable name");
            }
          } else fatal("invalid variable");
          
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
          if (!j.contains("variable")) fatal("missing nc variable");

          int ncid, varid, ncdims, nd;
          NC_SAFE_CALL( nc_open(filename0.c_str(), NC_NOWRITE, &ncid) );
          if (is_single_component()) {
            NC_SAFE_CALL( nc_inq_varid(ncid, j["variable"].template get<std::string>().c_str(), &varid) );
            j["varid"] = varid;
          } else {
            for (int i = 0; i < j["variable"].size(); i ++) {
              NC_SAFE_CALL( nc_inq_varid( ncid, j["variable"][i].template get<std::string>().c_str(), &varid) );
              j["varid"][i] = varid;
            }
          }
          NC_SAFE_CALL( nc_inq_varndims(ncid, varid, &ncdims) );

          // determin spatial dimensions
          int dimids[4]; 
          size_t dimlens[4];
          NC_SAFE_CALL( nc_inq_vardimid(ncid, varid, dimids) );
          for (int i = 0; i < ncdims; i ++) 
            NC_SAFE_CALL( nc_inq_dimlen(ncid, dimids[i], &dimlens[i]) );
          nc_close(ncid);

          if (ncdims == 4) { // 3 spatial dims + 1 time dimension; nd not required
            j["width"] = dimlens[3];
            j["height"] = dimlens[2];
            j["depth"] = dimlens[1];
            j["nd"] = 3;
          } else if (ncdims == 3) {
            if (j.contains("nd") && j["nd"].is_number()) { // nd required
              if (j["nd"] == 3) {
                j["width"] = dimlens[2];
                j["height"] = dimlens[1];
                j["depth"] = dimlens[0];
              } else if (j["nd"] == 2) {
                j["width"] = dimlens[2];
                j["height"] = dimlens[1];
              }
            } else fatal("missing/invalid nd");
          } else if (ncdims == 2) { // nd not required
            j["width"] = dimlens[1];
            j["height"] = dimlens[0];
            j["nd"] = 2;
          } else fatal("unsupported netcdf variable dimensionality");

          // determine number timesteps
          if (j.contains("n_timesteps") && j["n_timesteps"].is_number())
            j["n_timesteps"] = std::min(j["n_timesteps"].template get<size_t>(), j["filenames"].size());
          else 
            j["n_timesteps"] = j["filenames"].size();
#else
          fatal("FTK not compiled with NetCDF.");
#endif
        } else if (j["format"] == "h5") {
          fatal("array stream w/ h5 not implemented yet");
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
  std::vector<size_t> shape;
  if (j["nd"] == 2) {
    if (is_single_component()) 
      shape = std::vector<size_t>({j["width"], j["height"]});
    else 
      shape = std::vector<size_t>({n_components(), j["width"], j["height"]});
  } else {
    if (is_single_component()) 
      shape = std::vector<size_t>({j["width"], j["height"], j["depth"]});
    else 
      shape = std::vector<size_t>({n_components(), j["width"], j["height"], j["depth"]});
  }
  return shape;
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_file(int k)
{
  if (j["format"] == "float32") 
    return request_timestep_binary<float>(k);
  else if (j["format"] == "float64")
    return request_timestep_binary<double>(k);
  else if (j["format"] == "vti")
    return request_timestep_file_vti(k);
  else if (j["format"] == "netcdf")
    return request_timestep_file_nc(k);
  else if (j["format"] == "h5")
    return request_timestep_file_h5(k);
  else return ndarray<T>();
}

template <typename T>
template <typename T1>
ndarray<T> ndarray_stream<T>::request_timestep_binary(int k)
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
#if FTK_HAVE_VTK
  if (is_single_component())
    array.from_vtk_image_data_file(filename, j["variable"]);
  else {
    const int nv = n_components();
    std::vector<ftk::ndarray<T>> arrays(nv);
    for (int i = 0; i < nv; i ++)
      arrays[i].from_vtk_image_data_file(filename, j["variable"][i]);

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
    array.from_netcdf(filename, j["variable"]);
    array.reshape(shape()); // ncdims may not be equal to nd
  } else { // u, v, w in separate variables
    const int nv = n_components();
    std::vector<ftk::ndarray<T>> arrays(nv);
    for (int i = 0; i < nv; i ++)
      arrays[i].from_netcdf(filename, j["variable"][i]);

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
  fatal("h5 not yet implemented");
  return ndarray<T>();
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic(int k)
{
  if (j["name"] == "woven") 
    return request_timestep_synthetic_woven(k);
  else if (j["name"] == "double_gyre")
    return request_timestep_synthetic_double_gyre(k);
  else if (j["name"] == "merger")
    return request_timestep_synthetic_merger(k);
  else return ndarray<T>();
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_woven(int k) 
{
  const int nt = j["n_timesteps"];
  const T t = nt == 1 ? 0.0 : double(k)/(nt-1);
  return ftk::synthetic_woven_2D<T>(j["width"], j["height"], t);
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_merger(int k) 
{
  return ndarray<T>(); // TODO
}

template <typename T>
ndarray<T> ndarray_stream<T>::request_timestep_synthetic_double_gyre(int k) 
{
  return ndarray<T>(); // TODO
}

template <typename T>
void ndarray_stream<T>::start()
{
  if (!callback) 
    fatal("callback function not set");

  auto my_callback = callback;
  if (j.contains("temporal-smoothing-kernel")) {
    temporal_filter.set_gaussian_kernel(j["temporal-smoothing-kernel"], j["temporal-smoothing-kernel-size"]);
    temporal_filter.set_callback(callback);
    my_callback = [&](int k, ndarray<T>& array) {
      temporal_filter.push(array);
    };
  }

  for (int i = 0; i < j["n_timesteps"]; i ++) {
    ndarray<T> array;
    if (j["type"] == "synthetic") 
      array = request_timestep_synthetic(i);
    else if (j["type"] == "file")
      ndarray<T> array = request_timestep_file(i);
    my_callback(i, array);
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
