#ifndef __FTK_ERROR_HH
#define __FTK_ERROR_HH

#include <ftk/config.hh>
#include <execinfo.h>

namespace ftk {

enum {
  FTK_ERR_NOT_IMPLEMENTED = 1,
  FTK_ERR_FILE_NOT_FOUND = 1000,
  FTK_ERR_FILE_CANNOT_OPEN,
  FTK_ERR_FILE_CANNOT_WRITE,
  FTK_ERR_FILE_CANNOT_READ_EXPECTED_BYTES,
  FTK_ERR_FILE_UNRECOGNIZED_EXTENSION,
  FTK_ERR_FILE_FORMAT,
  FTK_ERR_FILE_FORMAT_AMIRA,
  FTK_ERR_NOT_BUILT_WITH_ADIOS2 = 2000,
  FTK_ERR_NOT_BUILT_WITH_ADIOS1,
  FTK_ERR_NOT_BUILT_WITH_BOOST,
  FTK_ERR_NOT_BUILT_WITH_CGAL,
  FTK_ERR_NOT_BUILT_WITH_CUDA,
  FTK_ERR_NOT_BUILT_WITH_GMP,
  FTK_ERR_NOT_BUILT_WITH_HDF5,
  FTK_ERR_NOT_BUILT_WITH_HIPSYCL,
  FTK_ERR_NOT_BUILT_WITH_SYCL,
  FTK_ERR_NOT_BUILT_WITH_KOKKOS,
  FTK_ERR_NOT_BUILT_WITH_LEVELDB,
  FTK_ERR_NOT_BUILT_WITH_METIS,
  FTK_ERR_NOT_BUILT_WITH_MPI,
  FTK_ERR_NOT_BUILT_WITH_MPSOLVE,
  FTK_ERR_NOT_BUILT_WITH_NETCDF,
  FTK_ERR_NOT_BUILT_WITH_OPENMP,
  FTK_ERR_NOT_BUILT_WITH_PARAVIEW,
  FTK_ERR_NOT_BUILT_WITH_PNETCDF,
  FTK_ERR_NOT_BUILT_WITH_PNG,
  FTK_ERR_NOT_BUILT_WITH_PYBIND11,
  FTK_ERR_NOT_BUILT_WITH_ROCKSDB,
  FTK_ERR_NOT_BUILT_WITH_QT5,
  FTK_ERR_NOT_BUILT_WITH_QT,
  FTK_ERR_NOT_BUILT_WITH_TBB,
  FTK_ERR_NOT_BUILT_WITH_VTK,
  FTK_ERR_NDARRAY_MULTIDIMENSIONAL_COMPONENTS = 3000, // only support one dim for components
  FTK_ERR_NDARRAY_UNSUPPORTED_DIMENSIONALITY,
  FTK_ERR_NDARRAY_RESHAPE_EMPTY,
  FTK_ERR_ACCELERATOR_UNSUPPORTED = 4000,
  FTK_ERR_THREAD_BACKEND_UNSUPPORTED = 5000,
  FTK_ERR_VTK_VARIABLE_NOT_FOUND = 6000,
  FTK_ERR_VTK_UNSUPPORTED_OUTPUT_FORMAT,
  FTK_ERR_NETCDF_MISSING_VARIABLE = 6500,
  FTK_ERR_ADIOS2 = 7000,
  FTK_ERR_ADIOS2_VARIABLE_NOT_FOUND,
  FTK_ERR_MESH_UNSUPPORTED_FORMAT = 8000,
  FTK_ERR_MESH_NONSIMPLICIAL, 
  FTK_ERR_MESH_EMPTY,
  FTK_ERR_UNKNOWN_OPTIONS = 10000
};

inline std::string err2str(int e)
{
  switch (e) {
  case FTK_ERR_NOT_IMPLEMENTED: return "not implemented yet";
  case FTK_ERR_FILE_NOT_FOUND: return "file not found";
  case FTK_ERR_FILE_CANNOT_OPEN: return "cannot open file";
  case FTK_ERR_FILE_CANNOT_WRITE: return "cannot write file";
  case FTK_ERR_FILE_CANNOT_READ_EXPECTED_BYTES: return "cannot read expected number of bytes";
  case FTK_ERR_FILE_FORMAT: return "file format error";
  case FTK_ERR_FILE_FORMAT_AMIRA: return "file format error with AmiraMesh data";
  case FTK_ERR_FILE_UNRECOGNIZED_EXTENSION: return "unrecognized file extension";
  case FTK_ERR_NOT_BUILT_WITH_ADIOS2: return "FTK not compiled with ADIOS2";
  case FTK_ERR_NOT_BUILT_WITH_ADIOS1: return "FTK not compiled with ADIOS1";
  case FTK_ERR_NOT_BUILT_WITH_BOOST: return "FTK not compiled with Boost";
  case FTK_ERR_NOT_BUILT_WITH_CGAL: return "FTK not compiled with CGAL";
  case FTK_ERR_NOT_BUILT_WITH_CUDA: return "FTK not compiled with CUDA";
  case FTK_ERR_NOT_BUILT_WITH_GMP: return "FTK not compiled with GMP";
  case FTK_ERR_NOT_BUILT_WITH_HDF5: return "FTK not compiled with HDF5";
  case FTK_ERR_NOT_BUILT_WITH_HIPSYCL: return "FTK not compiled with hipSYCL";
  case FTK_ERR_NOT_BUILT_WITH_KOKKOS: return "FTK not compiled with Kokkos";
  case FTK_ERR_NOT_BUILT_WITH_LEVELDB: return "FTK not compiled with LevelDB";
  case FTK_ERR_NOT_BUILT_WITH_METIS: return "FTK not compiled with Metis";
  case FTK_ERR_NOT_BUILT_WITH_MPI: return "FTK not compiled with MPI";
  case FTK_ERR_NOT_BUILT_WITH_MPSOLVE: return "FTK not compiled with MPSolve";
  case FTK_ERR_NOT_BUILT_WITH_NETCDF: return "FTK not compiled with NetCDF";
  case FTK_ERR_NOT_BUILT_WITH_OPENMP: return "FTK not compiled with OpenMP";
  case FTK_ERR_NOT_BUILT_WITH_PARAVIEW: return "FTK not compiled with ParaView";
  case FTK_ERR_NOT_BUILT_WITH_PNETCDF: return "FTK not compiled with Parallel-NetCDF";
  case FTK_ERR_NOT_BUILT_WITH_PNG: return "FTK not compiled with PNG";
  case FTK_ERR_NOT_BUILT_WITH_PYBIND11: return "FTK not compiled with PyBind11";
  case FTK_ERR_NOT_BUILT_WITH_ROCKSDB: return "FTK not compiled with RocksDB";
  case FTK_ERR_NOT_BUILT_WITH_QT5: return "FTK not compiled with Qt5";
  case FTK_ERR_NOT_BUILT_WITH_QT: return "FTK not compiled with Qt";
  case FTK_ERR_NOT_BUILT_WITH_TBB: return "FTK not compiled with TBB";
  case FTK_ERR_NOT_BUILT_WITH_VTK: return "FTK not compiled with VTK";
  case FTK_ERR_NDARRAY_MULTIDIMENSIONAL_COMPONENTS: return "FTK only supports one dim for components";
  case FTK_ERR_NDARRAY_UNSUPPORTED_DIMENSIONALITY: return "unsupported data dimensionality";
  case FTK_ERR_NDARRAY_RESHAPE_EMPTY: return "unable to reshape empty array";
  case FTK_ERR_ACCELERATOR_UNSUPPORTED: return "unsupported accelerator";
  case FTK_ERR_THREAD_BACKEND_UNSUPPORTED: return "unsupported thread backend";
  case FTK_ERR_VTK_VARIABLE_NOT_FOUND: return "VTK variable not found";
  case FTK_ERR_VTK_UNSUPPORTED_OUTPUT_FORMAT: return "unsupported vtk output format";
  case FTK_ERR_NETCDF_MISSING_VARIABLE: return "missing netcdf variable name(s)";
  case FTK_ERR_ADIOS2: return "adios2 error";
  case FTK_ERR_ADIOS2_VARIABLE_NOT_FOUND: return "adios2 variable not found";
  case FTK_ERR_MESH_UNSUPPORTED_FORMAT: return "unsupported mesh format";
  case FTK_ERR_MESH_NONSIMPLICIAL: return "unsupported nonsimplicial mesh";
  case FTK_ERR_MESH_EMPTY: return "empty mesh";
  default: return "unknown error: " + std::to_string(e);
  }
}

inline void print_backtrace()
{
  void *array[10];
  size_t size;
  char **strings;
  size_t i;

  size = backtrace (array, 10);
  strings = backtrace_symbols (array, size);

  printf ("Obtained %zd stack frames.\n", size);

  for (i = 0; i < size; i++)
    printf ("%s\n", strings[i]);

  free (strings);
}

inline void fatal(int err, std::string str = "")
{
  std::cerr << "[FTK FATAL] " << err2str(err);
  if (str.length()) std::cerr << ": " << str;
  std::cerr << std::endl;
  
  print_backtrace();
  exit(1);
}

inline void warn(int err, std::string str = "")
{
  std::cerr << "[FTK WARN] " << err2str(err);
  if (str.length()) std::cerr << ": " << str;
  std::cerr << std::endl;
}

inline void fatal(const std::string& str) {
  std::cerr << "[FTK FATAL] " << str << std::endl;
  
  print_backtrace();
  exit(1);
}

inline void warn(const std::string& str) {
  std::cerr << "[FTK WARN] " << str << std::endl;
}

}

#endif
