#ifndef _FTK_FILE_EXTENSIONS_HH
#define _FTK_FILE_EXTENSIONS_HH

#include <ftk/ftk_config.hh>
#include <ftk/utils/string.hh>

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

}

#endif
