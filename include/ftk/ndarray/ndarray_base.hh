#ifndef _FTK_NDARRAY_BASE_HH
#define _FTK_NDARRAY_BASE_HH

#include <ftk/config.hh>
#include <ftk/object.hh>
#include <ftk/error.hh>
#include <ftk/mesh/lattice.hh>
#include <ftk/io/util.hh>
#include <vector>
#include <array>
#include <numeric>
#include <tuple>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <random>

#if FTK_HAVE_VTK
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkDataReader.h>
#include <vtkNew.h>
#endif

#if FTK_HAVE_HDF5
#include <hdf5.h>
#endif

#if FTK_HAVE_ADIOS2
#include <adios2.h>
#endif

#if FTK_HAVE_ADIOS1
#include <adios.h>
#include <adios_read.h>
#include <adios_error.h>
#endif

namespace ftk {

enum {
  NDARRAY_TYPE_UNKNOWN = 0,
  NDARRAY_TYPE_FLOAT = 1,
  NDARRAY_TYPE_DOUBLE = 2,
  NDARRAY_TYPE_INT = 3
};

enum {
  NDARRAY_ADIOS2_STEPS_UNSPECIFIED = -1, 
  NDARRAY_ADIOS2_STEPS_ALL = -2
};

template <typename T> struct ndarray;

// the non-template base class for ndarray
struct ndarray_base {
  virtual int type() const = 0;

  virtual size_t size() const = 0;
  virtual bool empty() const = 0;

  size_t nd() const {return dims.size();}
  size_t dim(size_t i) const {return dims[i];}
  size_t shape(size_t i) const {return dim(i);}
  const std::vector<size_t> &shape() const {return dims;}
  size_t nelem() const { if (empty()) return 0; else return std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>()); }
  
  virtual void reshape(const std::vector<size_t> &dims_) = 0;
  void reshape(const std::vector<int>& dims);
  void reshape(size_t ndims, const size_t sizes[]);
  void reshape(const ndarray_base& array); //! copy shape from another array
  template <typename T> void reshape(const ndarray<T>& array); //! copy shape from another array
  
  size_t index(const std::vector<size_t>& idx) const;
  size_t index(const std::vector<int>& idx) const;
  
  template <typename uint=size_t>
  std::vector<uint> from_index(uint i) const {return lattice().from_integer(i);}

  lattice get_lattice() const;

public:
  virtual double get(size_t i0) const = 0;
  virtual double get(size_t i0, size_t i1) const = 0;
  virtual double get(size_t i0, size_t i1, size_t i2) const = 0;
  virtual double get(size_t i0, size_t i1, size_t i2, size_t i3) const = 0;
  virtual double get(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4) const = 0;
  virtual double get(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5) const = 0;
  virtual double get(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5, size_t i6) const = 0;

public:  
  void set_multicomponents(size_t c=1) {ncd = c;}
  size_t multicomponents() const {return ncd;}

  void set_has_time(bool b) { tv = b; }
  bool has_time() const { return tv; }

public: // binary i/o
  void read_binary_file(const std::string& filename);
  virtual void read_binary_file(FILE *fp) = 0;
  void to_binary_file(const std::string& filename);
  virtual void to_binary_file(FILE *fp) = 0;

public: // h5 i/o
  bool read_h5(const std::string& filename, const std::string& name);
#if FTK_HAVE_HDF5
  bool read_h5(hid_t fid, const std::string& name);
  virtual bool read_h5_did(hid_t did) = 0;
#endif

public: // adios2 i/o
  virtual void read_bp(
      const std::string& filename, 
      const std::string& varname, 
      int step = NDARRAY_ADIOS2_STEPS_UNSPECIFIED, 
      diy::mpi::communicator comm = MPI_COMM_WORLD) = 0;

#if FTK_HAVE_ADIOS2
  virtual void read_bp(
      adios2::IO &io, 
      adios2::Engine& reader, 
      const std::string &varname, 
      int step = NDARRAY_ADIOS2_STEPS_UNSPECIFIED) = 0; // read all
#endif

public: // adios1 io
  virtual bool read_bp_legacy(const std::string& filename, const std::string& varname, diy::mpi::communicator comm) = 0;

public: // vti i/o
  void read_vtk_image_data_file(const std::string& filename, const std::string array_name=std::string());
  virtual void read_vtk_image_data_file_sequence(const std::string& pattern) = 0;
#if FTK_HAVE_VTK
  virtual void from_vtk_image_data(vtkSmartPointer<vtkImageData> d, const std::string array_name=std::string()) = 0;
  virtual void from_vtu(vtkSmartPointer<vtkUnstructuredGrid> d, const std::string array_name=std::string()) = 0;
#endif

protected:
  std::vector<size_t> dims, s;
  size_t ncd = 0; // number of dimensions for components.  For 3D vector field, nd=4, ncd=1.  For 3D jacobian field, nd=5, ncd=2
  bool tv = false; // wheter the last dimension is time
};

////////
inline lattice ndarray_base::get_lattice() const {
  std::vector<size_t> st(nd(), 0), sz(dims);
  return lattice(st, sz);
}

inline void ndarray_base::reshape(const std::vector<int>& dims)
{
  std::vector<size_t> mydims;
  for (int i = 0; i < dims.size(); i ++)
    mydims.push_back(dims[i]);
  reshape(mydims);
}

inline void ndarray_base::reshape(size_t ndims, const size_t dims[])
{
  std::vector<size_t> mydims(dims, dims+ndims);
  reshape(mydims);
}

inline void ndarray_base::reshape(const ndarray_base& array)
{
  reshape(array.shape());
}

inline size_t ndarray_base::index(const std::vector<size_t>& idx) const {
  size_t i(idx[0]);
  for (size_t j = 1; j < nd(); j ++)
    i += idx[j] * s[j];
  return i;
}

inline size_t ndarray_base::index(const std::vector<int>& idx) const {
  std::vector<size_t> myidx(idx.begin(), idx.end());
  return index(myidx);
}

inline void ndarray_base::read_binary_file(const std::string& filename)
{
  FILE *fp = fopen(filename.c_str(), "rb");
  read_binary_file(fp);
  fclose(fp);
}

inline void ndarray_base::to_binary_file(const std::string& f)
{
  FILE *fp = fopen(f.c_str(), "wb");
  to_binary_file(fp);
  fclose(fp);
}

inline void ndarray_base::read_vtk_image_data_file(const std::string& filename, const std::string array_name)
{
#if FTK_HAVE_VTK
  vtkNew<vtkXMLImageDataReader> reader;
  reader->SetFileName(filename.c_str());
  reader->Update();
  from_vtk_image_data(reader->GetOutput(), array_name);
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
#endif
}

inline bool ndarray_base::read_h5(const std::string& filename, const std::string& name)
{
#if FTK_HAVE_HDF5
  auto fid = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fid < 0) return false; else {
    bool succ = read_h5(fid, name);
    H5Fclose(fid);
    return succ;
  }
#else 
  fatal(FTK_ERR_NOT_BUILT_WITH_HDF5);
  return false;
#endif
}

#if FTK_HAVE_HDF5
inline bool ndarray_base::read_h5(hid_t fid, const std::string& name)
{
  auto did = H5Dopen2(fid, name.c_str(), H5P_DEFAULT);
  if (did < 0) return false; else {
    bool succ = read_h5_did(did);
    H5Dclose(did);
    return succ;
  }
}
#endif

inline void ndarray_base::read_bp(const std::string& filename, const std::string& varname, int step, diy::mpi::communicator comm)
{
#if FTK_HAVE_ADIOS2
#if ADIOS2_USE_MPI
  adios2::ADIOS adios(comm);
#else
  adios2::ADIOS adios;
#endif
  adios2::IO io = adios.DeclareIO("BPReader");
  adios2::Engine reader = io.Open(filename, adios2::Mode::Read); // , MPI_COMM_SELF);
  
  read_bp(io, reader, varname, step);
  reader.Close();
  
  // empty array; try legacy reader
  if (empty()) {
#if FTK_HAVE_ADIOS1
    read_bp_legacy(filename, varname, comm);
#else
    throw FTK_ERR_ADIOS2;
#endif
  }
  
  // if (empty()) read_bp_legacy(filename, varname, comm); 
#else
  warn(FTK_ERR_NOT_BUILT_WITH_ADIOS2);
  read_bp_legacy(filename, varname, comm);
#endif
}

} // namespace ftk

#endif
