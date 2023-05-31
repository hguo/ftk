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

#if FTK_HAVE_NETCDF
#include <netcdf.h>
#include <netcdf_meta.h>
#if NC_HAS_PARALLEL
#include <netcdf_par.h>
#endif
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
  virtual size_t elem_size() const = 0;
 
  virtual const void* pdata() const = 0;
  virtual void* pdata() = 0;

  virtual void reshape(const std::vector<size_t> &dims_) = 0;
  void reshape(const std::vector<int>& dims);
  void reshape(size_t ndims, const size_t sizes[]);
  void reshape(const ndarray_base& array); //! copy shape from another array
  template <typename T> void reshape(const ndarray<T>& array); //! copy shape from another array
  
  size_t index(const std::vector<size_t>& idx) const;
  size_t index(const std::vector<int>& idx) const;
  size_t index(const size_t idx[]) const;
  
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
  void make_multicomponents(); // make a non-multicomponent array to an array with 1 component
  size_t multicomponents() const {return ncd;}
  size_t ncomponents() const;

  void set_has_time(bool b) { tv = b; }
  bool has_time() const { return tv; }

public: // binary i/o
  void read_binary_file(const std::string& filename);
  virtual void read_binary_file(FILE *fp) = 0;
  void to_binary_file(const std::string& filename);
  virtual void to_binary_file(FILE *fp) = 0;

public: // netcdf
  void read_netcdf(const std::string& filename, const std::string& varname, const size_t starts[], const size_t sizes[], diy::mpi::communicator comm=MPI_COMM_WORLD);
  void read_netcdf(int ncid, const std::string& varname, const size_t starts[], const size_t sizes[], diy::mpi::communicator comm=MPI_COMM_WORLD);
  void read_netcdf(int ncid, int varid, int ndims, const size_t starts[], const size_t sizes[], diy::mpi::communicator comm=MPI_COMM_WORLD);
  void read_netcdf(int ncid, int varid, const size_t starts[], const size_t sizes[], diy::mpi::communicator comm=MPI_COMM_WORLD);
  void read_netcdf(const std::string& filename, const std::string& varname, diy::mpi::communicator comm=MPI_COMM_WORLD);
  void read_netcdf(int ncid, const std::string& varname, diy::mpi::communicator comm=MPI_COMM_WORLD);
  void read_netcdf(int ncid, int varid, diy::mpi::communicator comm=MPI_COMM_WORLD);
  void read_netcdf_slice(const std::string& filename, const std::string& varname, int k, diy::mpi::communicator comm=MPI_COMM_WORLD);
  // void to_netcdf(int ncid, const std::string& varname);
  // void to_netcdf(int ncid, int varid);
  void to_netcdf(int ncid, int varid, const size_t starts[], const size_t sizes[]) const;
  void to_netcdf(int ncid, int varid) const;
  void to_netcdf_multivariate(int ncid, int varids[]) const;
  void to_netcdf_unlimited_time(int ncid, int varid) const;
  void to_netcdf_multivariate_unlimited_time(int ncid, int varids[]) const;
  
  virtual int nc_datatype() const = 0;

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

public: // vtk data array
#if FTK_HAVE_VTK
  vtkSmartPointer<vtkDataArray> to_vtk_data_array(std::string varname=std::string()) const; 
  virtual int vtk_data_type() const = 0;
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

inline size_t ndarray_base::ncomponents() const {
  size_t rtn = 1;
  for (size_t i = 0; i < multicomponents(); i ++)
    rtn *= dims[i];
  return rtn;
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

inline size_t ndarray_base::index(const size_t idx[]) const {
  size_t i(idx[0]);
  for (size_t j = 1; j < nd(); j ++)
    i += idx[j] * s[j];
  return i;
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

inline void ndarray_base::make_multicomponents()
{
  std::vector<size_t> s = shape();
  s.insert(s.begin(), 1);
  reshape(s);
  set_multicomponents();
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

#if FTK_HAVE_VTK
inline vtkSmartPointer<vtkDataArray> ndarray_base::to_vtk_data_array(std::string varname) const
{
  vtkSmartPointer<vtkDataArray> d = vtkDataArray::CreateDataArray(this->vtk_data_type());
  if (varname.length() > 0)
    d->SetName( varname.c_str() );

  // fprintf(stderr, "to_vtk_data_array, ncd=%zu\n", ncd);

  if (ncd == 1) {
    d->SetNumberOfComponents(shape(0));
    d->SetNumberOfTuples( std::accumulate(dims.begin()+1, dims.end(), 1, std::multiplies<size_t>()) );
  }
  else if (ncd == 0) {
    d->SetNumberOfComponents(1);
    d->SetNumberOfTuples(nelem());
  } else {
    fatal(FTK_ERR_NDARRAY_MULTIDIMENSIONAL_COMPONENTS);
  }
  memcpy(d->GetVoidPointer(0), this->pdata(), elem_size() * nelem()); // nelem());
  return d;
}
#endif

inline void ndarray_base::read_netcdf(const std::string& filename, const std::string& varname, diy::mpi::communicator comm)
{
#if FTK_HAVE_NETCDF
  int ncid, varid;
#if NC_HAS_PARALLEL
  int rtn = nc_open_par(filename.c_str(), NC_NOWRITE, comm, MPI_INFO_NULL, &ncid);
  if (rtn != NC_NOERR)
    NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
#else
  NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
#endif
  NC_SAFE_CALL( nc_inq_varid(ncid, varname.c_str(), &varid) );
  read_netcdf(ncid, varid, comm);
  NC_SAFE_CALL( nc_close(ncid) );
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

inline void ndarray_base::read_netcdf(int ncid, int varid, int ndims, const size_t starts[], const size_t sizes[], diy::mpi::communicator comm)
{
#if FTK_HAVE_NETCDF
  std::vector<size_t> mysizes(sizes, sizes+ndims);
  std::reverse(mysizes.begin(), mysizes.end());
  reshape(mysizes);

  if (nc_datatype() == NC_INT) {
    NC_SAFE_CALL( nc_get_vara_int(ncid, varid, starts, sizes, (int*)pdata()) );
  } else if (nc_datatype() == NC_FLOAT) {
    NC_SAFE_CALL( nc_get_vara_float(ncid, varid, starts, sizes, (float*)pdata()) );
  } else if (nc_datatype() == NC_DOUBLE) {
    NC_SAFE_CALL( nc_get_vara_double(ncid, varid, starts, sizes, (double*)pdata()) );
  } else if (nc_datatype() == NC_UINT) {
    NC_SAFE_CALL( nc_get_vara_uint(ncid, varid, starts, sizes, (unsigned int*)pdata()) );
  } else 
    fatal(FTK_ERR_NOT_IMPLEMENTED);
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

inline void ndarray_base::to_netcdf(int ncid, int varid) const
{
  std::vector<size_t> starts(dims.size(), 0), sizes(dims);
  std::reverse(sizes.begin(), sizes.end());

  to_netcdf(ncid, varid, &starts[0], &sizes[0]);
}

inline void ndarray_base::to_netcdf(int ncid, int varid, const size_t st[], const size_t sz[]) const
{
#ifdef FTK_HAVE_NETCDF
  fprintf(stderr, "st=%zu, %zu, %zu, %zu, sz=%zu, %zu, %zu, %zu\n", 
      st[0], st[1], st[2], st[3], sz[0], sz[1], sz[2], sz[3]);
  
  if (nc_datatype() == NC_DOUBLE) {
    NC_SAFE_CALL( nc_put_vara_double(ncid, varid, st, sz, (double*)pdata()) );
  } else if (nc_datatype() == NC_FLOAT) {
    NC_SAFE_CALL( nc_put_vara_float(ncid, varid, st, sz, (float*)pdata()) );
  } else if (nc_datatype() == NC_INT) {
    NC_SAFE_CALL( nc_put_vara_int(ncid, varid, st, sz, (int*)pdata()) );
  } else 
    fatal(FTK_ERR_NOT_IMPLEMENTED);
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

#if 0
inline void ndarray_base::to_netcdf_multivariate(int ncid, int varids[]) const
{
  const size_t nv = dims[0], ndims = nd()-1;
  std::vector<size_t> d(dims.begin()+1, dims.end());

  for (int i = 0; i < nv; i ++) {
    ndarray<T> subarray(d);
    for (size_t j = 0; j < subarray.nelem(); j ++) 
      subarray[j] = p[j*nv + i];
    subarray.to_netcdf(ncid, varids[i]);
  }
}
#endif

inline void ndarray_base::to_netcdf_unlimited_time(int ncid, int varid) const
{
  std::vector<size_t> starts(dims.size()+1, 0), sizes(dims);
  sizes.push_back(1);
  std::reverse(sizes.begin(), sizes.end());
 
  // fprintf(stderr, "starts={%zu, %zu, %zu}, sizes={%zu, %zu, %zu}\n", 
  //     starts[0], starts[1], starts[2], sizes[0], sizes[1], sizes[2]);

  to_netcdf(ncid, varid, &starts[0], &sizes[0]);
}

#if 0
inline void ndarray_base::to_netcdf_multivariate_unlimited_time(int ncid, int varids[]) const
{
  const size_t nv = dims[0], ndims = nd()-1;
  std::vector<size_t> d(dims.begin()+1, dims.end());

  for (int i = 0; i < nv; i ++) {
    ndarray<T> subarray(d);
    for (size_t j = 0; j < subarray.nelem(); j ++) 
      subarray[j] = p[j*nv + i];
    subarray.to_netcdf_unlimited_time(ncid, varids[i]);
  }
}
#endif

inline void ndarray_base::read_netcdf(int ncid, int varid, const size_t starts[], const size_t sizes[], diy::mpi::communicator comm)
{
#ifdef FTK_HAVE_NETCDF
  int ndims;
  NC_SAFE_CALL( nc_inq_varndims(ncid, varid, &ndims) );

  std::vector<size_t> mysizes(sizes, sizes+ndims);
  std::reverse(mysizes.begin(), mysizes.end());
  reshape(mysizes);

  read_netcdf(ncid, varid, ndims, starts, sizes, comm);
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

inline void ndarray_base::read_netcdf(int ncid, int varid, diy::mpi::communicator comm)
{
#ifdef FTK_HAVE_NETCDF
  int ndims;
  int dimids[4];
  size_t starts[4] = {0}, sizes[4] = {0};

  NC_SAFE_CALL( nc_inq_varndims(ncid, varid, &ndims) );
  NC_SAFE_CALL( nc_inq_vardimid(ncid, varid, dimids) );

  for (int i = 0; i < ndims; i ++)
    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimids[i], &sizes[i]) );
  
  read_netcdf(ncid, varid, starts, sizes, comm);
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

inline void ndarray_base::read_netcdf(int ncid, const std::string& varname, diy::mpi::communicator comm)
{
#ifdef FTK_HAVE_NETCDF
  int varid;
  NC_SAFE_CALL( nc_inq_varid(ncid, varname.c_str(), &varid) );
  read_netcdf(ncid, varid, comm);
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

inline void ndarray_base::read_netcdf(int ncid, const std::string& varname, const size_t starts[], const size_t sizes[], diy::mpi::communicator comm)
{
#ifdef FTK_HAVE_NETCDF
  int varid;
  NC_SAFE_CALL( nc_inq_varid(ncid, varname.c_str(), &varid) );
  read_netcdf(ncid, varid, starts, sizes, comm);
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

inline void ndarray_base::read_netcdf(const std::string& filename, const std::string& varname, const size_t starts[], const size_t sizes[], diy::mpi::communicator comm)
{
#ifdef FTK_HAVE_NETCDF
  int ncid, varid;
#if NC_HAS_PARALLEL
  int rtn = nc_open_par(filename.c_str(), NC_NOWRITE, comm, MPI_INFO_NULL, &ncid);
  if (rtn != NC_NOERR)
    NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
#else
  NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
#endif

  NC_SAFE_CALL( nc_inq_varid(ncid, varname.c_str(), &varid) );
  read_netcdf(ncid, varid, starts, sizes, comm);
  NC_SAFE_CALL( nc_close(ncid) );
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

} // namespace ftk

#endif
