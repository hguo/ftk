#ifndef _FTK_NDARRAY_STREAM_WRITER_HH
#define _FTK_NDARRAY_STREAM_WRITER_HH

#include <ftk/object.hh>
#include <ftk/ndarray/stream.hh>

namespace ftk {
using nlohmann::json;

template <typename T=double>
struct ndarray_writer : public object {
  void configure(const json& j_);
  // JSON specificications:
  // required fields: 
  //  - nd, number. dimensionality of the data
  //  - filename, string, e.g. "tornado-%04d.nc".  Should be able to expand to file names with sprintf(3)
  //  - format (nc|vti|float32|float64)
  // format-specific fields: 
  //  - nc (NetCDF)
  //    - variable (required), string or array of strings.  Possible to use the same
  //      specification in the input stream config.  Will be converted to array of 
  //      strings internally.  Will write one or multiple variables.
  //    - dimensions (optional), array of strings; by default use ["x", "y"] and 
  //      ["x", "y", "z"] for 2D and 3D data, respectively
  //    - unlimited_time (optional, by default true), bool.  If true, include time 
  //      in dimensions, e.g. (time,z,y,x), in netcdf.  
  //  - vti (vtkImageData w/ vtkXMLImageDataWriter)
  //    - variable (required), string or array of strings.  Possible to use the same
  //      specification in the input stream config.  
  //      - string: write one single variable, assuming the data is single-component
  //      - array of strings: write one single variable, assuming the data is multi-
  //        component.  The size of the array should be the exact number of components
  //    - separate_variables (optional, by default false), bool
  //      - write separate variables instead of one single variable
  
  void consume(ndarray_stream<T>&);

protected:
  void write(int k, const ndarray<T> &data);
  void write_netcdf(int k, const ndarray<T> &data);
  void write_vti(int k, const ndarray<T> &data);
  std::string filename(int k) const;

public:
  static std::string filename(const std::string& pattern, int k);

private:
  json j;
};

//////////////////
template <typename T>
void ndarray_writer<T>::configure(const json& j_)
{
  j = j_; 

  // sanity check
  if (j.contains("nd")) {
    if (j["nd"].is_number()) {
      // OK
    } else fatal("invalid nd");
  } else fatal("missing nd");

  if (j.contains("filename")) {
    if (j["filename"].is_string()) {
      // TODO: check if the pattern can be expanded with sprintf
    } else fatal("invalid filename");
  } else fatal ("missing filename");

  if (j.contains("format")) {
    if (j["format"].is_string()) {
      const std::string format = j["format"];
      const std::set<std::string> valid = {"nc", "vti", "float32", "float64"};
      if (valid.find(format) == valid.end())
        fatal("invalid format");
    } else fatal("invalid foramt");
  } else fatal("missing format");

  if (j.contains("n_components")) {
    if (j["n_components"].is_number()) {
      // OK
    } else fatal("invalid n_components");
  } else j["n_components"] = 1;

  if (j.contains("dimensions")) {
    if (j["dimensions"].is_array()) {
      auto jd = j["dimensions"];
      if (jd.size() < 2 || jd.size() > 3) 
        fatal("invalid number of dimensions");
      for (int i = 0; i < jd.size(); i ++) {
        if (jd[i].is_string()) {
          // OK
        } else fatal("invalid dimensions");
      }
    } else fatal("invalid dimensions");
  } else {
    if (j["nd"] == 2) j["dimensions"] = {"x", "y"};
    else j["dimensions"] = {"x", "y", "z"};
  }

  if (j.contains("unlimited_time")) {
    if (j["format"] != "nc") warn("ignoring unlimited_time");
    else {
      if (j["unlimited_time"].is_boolean()) {
        // OK
      } else 
        fatal("invalid unlimited_time");
    }
  } else {
    if (j["format"] == "nc") j["unlimited_time"] = true;
  }

  if (j.contains("variable")) {
    if (j["format"] == "float32" || j["format"] == "float64") 
      warn("variable ignored");

    if (j["variable"].is_string()) {
      // convert to an array of one single string
      std::vector<std::string> vars;
      vars.push_back(j["variable"]);
      j["variable"] = vars;
    } else if (j["variable"].is_array()) {
      for (int i = 0; i < j["variable"].size(); i ++) 
        if (j["variable"][i].is_string()) {
          // OK
        } else 
          fatal("invalid variable");
    } else 
      fatal("invalid variable");
  } else {
    if (j["foramt"] == "nc") 
      fatal("variable name missing");
  }

  if (j.contains("separate_variables")) {
    if (j["format"] == "float32" || j["format"] == "float64") 
      warn("separate_variables ignored");

    if (j["separate_variables"].is_boolean()) {
      if (j["variable"].is_string()) 
        warn("separate_variables ignored");
    } else 
      fatal("invalid separate_variables");
  }
}

template <typename T>
std::string ndarray_writer<T>::filename(int k) const {
  const std::string pattern = j["filename"];
  return filename(pattern, k);
}

template <typename T>
std::string ndarray_writer<T>::filename(const std::string& pattern, int k)
{
  ssize_t size = snprintf(NULL, 0, pattern.c_str(), k);
  char *buf = (char*)malloc(size + 1);
  snprintf(buf, size + 1, pattern.c_str(), k);
  const std::string filename(buf);
  free(buf);
  return filename;
}

template <typename T>
void ndarray_writer<T>::write(int k, const ndarray<T> &data)
{
  const std::string f = this->filename(k);
  const std::string format = j["format"];
  if (format == "float32") data.template to_binary_file2<float>(f);
  else if (format == "float64") data.template to_binary_file2<double>(f);
  else if (format == "nc") write_netcdf(k, data);
  else if (format == "vti") write_vti(k, data);
}

template <typename T>
void ndarray_writer<T>::write_vti(int k, const ndarray<T>& data)
{
  const std::string filename = this->filename(k);
  const int nv = j["variable"].size();
  const bool multicomponent = nv > 1;
  data.to_vtk_image_data_file(filename, multicomponent);
}

template <typename T>
void ndarray_writer<T>::write_netcdf(int k, const ndarray<T> &data)
{
#if FTK_HAVE_NETCDF
  const std::string filename = this->filename(k);

  const int nd = j["dimensions"].size();
  const int nv = j["variable"].size();

  int ncid, varids[nv], dimids[4]; // max dim is 4
  int ncndims = 0; // temp variable to index new dim ids
  const bool unlimited_time = j["unlimited_time"];
  
  NC_SAFE_CALL( nc_create(filename.c_str(), NC_CLOBBER | NC_64BIT_OFFSET, &ncid) );

  // dimensions
  if (unlimited_time) 
    NC_SAFE_CALL( nc_def_dim(ncid, "time", NC_UNLIMITED, &dimids[ncndims++]) );

  for (int i = 0; i < nd; i ++) {
    const std::string dimname = j["dimensions"][nd-i-1];
    NC_SAFE_CALL( nc_def_dim(ncid, dimname.c_str(), data.dim(nd-i-1), &dimids[ncndims++]) );
  }

  // variable(s)
  for (int i = 0; i < nv; i ++) {
    const std::string var_name = j["variable"][i];
    NC_SAFE_CALL( nc_def_var(ncid, var_name.c_str(), data.nc_datatype(), ncndims, dimids, &varids[i]) );
  }
  
  nc_enddef(ncid);

  // write data
  for (int i = 0; i < nv; i ++)
    if (unlimited_time)
      data.to_netcdf_unlimited_time(ncid, varids[i]);
    else
      data.to_netcdf(ncid, varids[i]);

  NC_SAFE_CALL( nc_close(ncid) );
#else
  fatal("FTK not compiled with netcdf");
#endif
}

template <typename T>
void ndarray_writer<T>::consume(ndarray_stream<T>& stream)
{
  stream.set_callback([&](int k, ndarray<T> data) {
    write(k, data);
  });
}

}

#endif
