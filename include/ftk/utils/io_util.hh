#ifndef _FTK_IO_UTIL_HH
#define _FTK_IO_UTIL_HH

#include <ftk/config.hh>
#include <ftk/utils/string.hh>
#include <iostream>

namespace ftk {

enum {
  FILE_EXT_NULL = 0,
  FILE_EXT_BIN = 1,
  FILE_EXT_JSON,
  FILE_EXT_NETCDF,
  FILE_EXT_HDF5,
  FILE_EXT_BP,  // adios2
  FILE_EXT_NUMPY, // numpy
  FILE_EXT_PNG,
  FILE_EXT_VTI, // vtk xml image data
  FILE_EXT_VTP, // vtk xml poly data
  FILE_EXT_VTU, // vtk xml unstructured grid data
  FILE_EXT_VTK, // legacy vtk format
  FILE_EXT_PLY, // surface
  FILE_EXT_STL  // surface
};

static inline std::string series_filename(
    const std::string& pattern, int k)
{
  ssize_t size = snprintf(NULL, 0, pattern.c_str(), k);
  char *buf = (char*)malloc(size + 1);
  snprintf(buf, size + 1, pattern.c_str(), k);
  const std::string filename(buf);
  free(buf);
  return filename;
}

static bool file_exists(const std::string& filename) {
  std::ifstream f(filename);
  return f.good();
}

static bool file_not_exists(const std::string& filename) { 
  return !file_exists(filename); 
}

static inline int file_extension(const std::string& f)
{
  auto m = [f](std::string e) { return ends_with_lower(f, e); };

  if (m("bin") || m("binary"))
    return FILE_EXT_BIN;
  else if (m("nc") || m("netcdf"))
    return FILE_EXT_NETCDF;
  else if (m("h5") || m("hdf5"))
    return FILE_EXT_HDF5;
  else if (m("bp") || m("adios2"))
    return FILE_EXT_HDF5;
  else if (m("npy") || m("numpy"))
    return FILE_EXT_NUMPY;
  else if (m("png"))
    return FILE_EXT_PNG;
  else if (m("vti"))
    return FILE_EXT_VTI;
  else if (m("vtp"))
    return FILE_EXT_VTP;
  else if (m("vtu"))
    return FILE_EXT_VTU;
  else if (m("vtk"))
    return FILE_EXT_VTK;
  else if (m("ply"))
    return FILE_EXT_PLY;
  else if (m("stl"))
    return FILE_EXT_STL;
  else 
    return FILE_EXT_NULL;
}

static inline int file_extension(const std::string& filename, const std::string& format)
{
  if (format == "auto" || format.empty())
    return file_extension(filename);
  else 
    return file_extension(format);
}

} // namespace ftk

#endif
