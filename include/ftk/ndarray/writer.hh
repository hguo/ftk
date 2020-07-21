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
  //      specification in the input stream config.
  //      - string: write one single variable, assuming the data is single-component
  //      - array of strings: write one single variable, assuming the data is multi-
  //        component.  The size of the array should be the exact number of components
  //    - dimension_names (required), e.g. use ["x", "y"] and ["x", "y", "z"] for 2D 
  //      and 3D data, respectively
  //    - time (optional, by default false), bool, number, or string
  //      - bool: if yes, include time in timension_names, e.g. (time,z,y,x), in 
  //        netcdf.  Use "time" as the dimension name, and use timestep as the time value
  //      - number: include time in dimension_names in netcdf.  Use "time" as the 
  //        dimension name, and use the designated value as the time value
  //      - string: include time and use the designated name as the dimension name
  //    - separate_variables (optional, by default false), bool
  //      - write separate variables instead of one single variable
  //  - vti (vtkImageData w/ vtkXMLImageDataWriter)
  //    - variable (required), see "variable" in NetCDF section above
  //    - separate_variables, see "separate_variables" in NetCDF section above
  
  void consume(ndarray_stream<T>&);

protected:
  void write(int k, const ndarray<T> &data);
  void write_netcdf(int k, const ndarray<T> &data);
  std::string filename(int k) const;

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

  if (j.contains("dimension_names")) {
    if (j["dimension_names"].is_array()) {
      auto jd = j["dimension_names"];
      if (jd.size() < 2 || jd.size() > 3) 
        fatal("invalid number of dimension_names");
      for (int i = 0; i < jd.size(); i ++) {
        if (jd[i].is_string()) {
          // OK
        } else fatal("invalid dimension_names");
      }
    } else fatal("invalid dimension_names");
  } else {
    if (j["nd"] == 2) j["dimension_names"] = {"x", "y"};
    else j["dimension_names"] = {"x", "y", "z"};
  }

  if (j.contains("variable")) {
    if (j["format"] == "float32" || j["format"] == "float64") 
      warn("variable ignored");

    if (j["variable"].is_string()) {
      // OK
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
}

template <typename T>
void ndarray_writer<T>::write_netcdf(int k, const ndarray<T> &data)
{
#if FTK_HAVE_NETCDF
  const std::string filename = this->filename(k);

  const int nd = j["dimension_names"];
  int nv;
  if (j["variable"].is_string()) nv = 1;
  else if (j["variable"].is_array()) nv = j["variable"].size();
  else fatal("missing netcdf variables");

  int ncid, varids[nv], dimids[nd];

  NC_SAFE_CALL( nc_create(filename.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid) );
  for (int i = 0; i < nd; i ++) {
    const std::string dimname = j["dimension_names"][i];
    NC_SAFE_CALL( nc_def_dim(ncid, dimname.c_str(), size, dimids[i]) );
  }

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
