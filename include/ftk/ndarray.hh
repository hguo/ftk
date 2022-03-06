#ifndef _FTK_NDARRAY_HH
#define _FTK_NDARRAY_HH

#include <ftk/config.hh>
#include <ftk/ndarray/ndarray_base.hh>
#include <ftk/numeric/clamp.hh>
#include <ftk/basic/murmurhash2.hh>

#if FTK_HAVE_CUDA
// #include <cuda.h>
// #include <cuda_runtime.h>
#endif

#if FTK_HAVE_MPI
#include <mpi.h>
#include <ftk/external/bil/bil.h>
#endif

#if FTK_HAVE_NETCDF
#include <netcdf.h>
#include <netcdf_meta.h>
#if NC_HAS_PARALLEL
#include <netcdf_par.h>
#endif
#endif

#if FTK_HAVE_PNETCDF
#include <pnetcdf.h>
#endif

#if FTK_HAVE_PYBIND11
#include <pybind11/numpy.h>
#endif

namespace ftk {

template <typename T>
struct ndarray : public ndarray_base {
  friend class diy::Serialization<ndarray<T>>;

  int type() const;

  unsigned int hash() const;

  ndarray() {}
  ndarray(const std::vector<size_t> &dims) {reshape(dims);}
  ndarray(const std::vector<size_t> &dims, T val) {reshape(dims, val);}
  ndarray(const lattice& l) {reshape(l.sizes());}
  ndarray(const T *a, const std::vector<size_t> &shape);
  
  template <typename T1> ndarray(const ndarray<T1>& array1) { from_array<T1>(array1); }
  ndarray(const ndarray<T>& a) { dims = a.dims; s = a.s; ncd = a.ncd; tv = a.tv; p = a.p; }
  
  template <typename T1> ndarray<T>& operator=(const ndarray<T1>& array1) { from_array<T1>(array1); return *this; }
  ndarray<T>& operator=(const ndarray<T>& a) { dims = a.dims; s = a.s; ncd = a.ncd; tv = a.tv; p = a.p; return *this; }

  std::ostream& print(std::ostream& os) const;
  std::ostream& print_shape(std::ostream& os) const;

  size_t size() const {return p.size();}
  bool empty() const  {return p.empty();}

  void fill(T value); //! fill with a constant value
  void fill(const std::vector<T>& values); //! fill values with std::vector
  void fill(const std::vector<std::vector<T>>& values); //! fill values

  const std::vector<T>& std_vector() const {return p;}

  const T* data() const {return p.data();}
  T* data() {return p.data();}

  void swap(ndarray& x);

  void reshape(const std::vector<size_t> &dims_);
  void reshape(const std::vector<size_t> &dims, T val);
  template <typename I> void reshape(const int ndims, const I sz[]);
  template <typename T1> void reshape(const ndarray<T1>& array); //! copy shape from another array

  void reshape(size_t n0) {reshape(std::vector<size_t>({n0}));}
  void reshape(size_t n0, size_t n1) {reshape({n0, n1});}
  void reshape(size_t n0, size_t n1, size_t n2) {reshape({n0, n1, n2});}
  void reshape(size_t n0, size_t n1, size_t n2, size_t n3) {reshape({n0, n1, n2, n3});}
  void reshape(size_t n0, size_t n1, size_t n2, size_t n3, size_t n4) {reshape({n0, n1, n2, n3, n4});}
  void reshape(size_t n0, size_t n1, size_t n2, size_t n3, size_t n4, size_t n5) {reshape({n0, n1, n2, n3, n4, n5});}
  void reshape(size_t n0, size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6) {reshape({n0, n1, n2, n3, n4, n5, n6});}

  ndarray<T> slice(const lattice&) const;
  ndarray<T> slice(const std::vector<size_t>& starts, const std::vector<size_t> &sizes) const;

  ndarray<T> slice_time(size_t t) const; // assuming the last dimension is time, return an (n-1)-dimensional slice
  std::vector<ndarray<T>> slice_time() const; // slice_time for all timesteps

  // merge multiple arrays into a multicomponent array
  static ndarray<T> concat(const std::vector<ndarray<T>>& arrays);
  static ndarray<T> stack(const std::vector<ndarray<T>>& arrays);

  T& at(const std::vector<size_t>& idx) {return p[index(idx)];}
  const T& at(const std::vector<size_t>& idx) const {return p[index(idx)];}
  
  T& at(const std::vector<int>& idx) {return p[index(idx)];}
  const T& at(const std::vector<int>& idx) const {return p[index(idx)];}

  T& at(size_t i0) {return p[i0];}
  T& at(size_t i0, size_t i1) {return p[i0+i1*s[1]];}
  T& at(size_t i0, size_t i1, size_t i2) {return p[i0+i1*s[1]+i2*s[2]];}
  T& at(size_t i0, size_t i1, size_t i2, size_t i3) {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]];}
  T& at(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4) {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]+i4*s[4]];}
  T& at(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5) {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]+i4*s[4]+i5*s[5]];}
  T& at(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5, size_t i6) {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]+i4*s[4]+i5*s[5]+i6*s[6]];}

  const T& at(size_t i0) const {return p[i0];}
  const T& at(size_t i0, size_t i1) const {return p[i0+i1*s[1]];}
  const T& at(size_t i0, size_t i1, size_t i2) const {return p[i0+i1*s[1]+i2*s[2]];}
  const T& at(size_t i0, size_t i1, size_t i2, size_t i3) const {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]];}
  const T& at(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4) const {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]+i4*s[4]];}
  const T& at(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5) const {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]+i4*s[4]+i5*s[5]];}
  const T& at(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5, size_t i6) const {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]+i4*s[4]+i5*s[5]+i6*s[6]];}

  double get(size_t i0) const { return at(i0); }
  double get(size_t i0, size_t i1) const { return at(i0, i1); }
  double get(size_t i0, size_t i1, size_t i2) const { return at(i0, i1, i2); }
  double get(size_t i0, size_t i1, size_t i2, size_t i3) const { return at(i0, i1, i2, i3); }
  double get(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4) const { return at(i0, i1, i2, i3, i4); }
  double get(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5) const { return at(i0, i1, i2, i3, i4, i5); }
  double get(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5, size_t i6) const { return at(i0, i1, i2, i3, i4, i5, i6); }

  T& operator()(const std::vector<size_t>& idx) {return at(idx);}
  T& operator()(const std::vector<int>& idx) {return at(idx);}
  T& operator()(size_t i0) {return p[i0];}
  T& operator()(size_t i0, size_t i1) {return p[i0+i1*s[1]];}
  T& operator()(size_t i0, size_t i1, size_t i2) {return p[i0+i1*s[1]+i2*s[2]];}
  T& operator()(size_t i0, size_t i1, size_t i2, size_t i3) {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]];}
  T& operator()(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4) {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]+i4*s[4]];}
  T& operator()(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5) {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]+i4*s[4]+i5*s[5]];}
  T& operator()(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5, size_t i6) {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]+i4*s[4]+i5*s[5]+i6*s[6]];}
  
  T& operator()(const std::vector<size_t>& idx) const {return at(idx);}
  T& operator()(const std::vector<int>& idx) const {return at(idx);}
  const T& operator()(size_t i0) const {return p[i0];}
  const T& operator()(size_t i0, size_t i1) const {return p[i0+i1*s[1]];}
  const T& operator()(size_t i0, size_t i1, size_t i2) const {return p[i0+i1*s[1]+i2*s[2]];}
  const T& operator()(size_t i0, size_t i1, size_t i2, size_t i3) const {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]];}
  const T& operator()(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4) const {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]+i4*s[4]];}
  const T& operator()(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5) const {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]+i4*s[4]+i5*s[5]];}
  const T& operator()(size_t i0, size_t i1, size_t i2, size_t i3, size_t i4, size_t i5, size_t i6) const {return p[i0+i1*s[1]+i2*s[2]+i3*s[3]+i4*s[4]+i5*s[5]+i6*s[6]];}

  friend std::ostream& operator<<(std::ostream& os, const ndarray<T>& arr) {arr.print(os); return os;}
  friend bool operator==(const ndarray<T>& lhs, const ndarray<T>& rhs) {return lhs.dims == rhs.dims && lhs.p == rhs.p;}

  ndarray<T>& operator+=(const ndarray<T>& x);
  ndarray<T>& operator-=(const ndarray<T>& x);
  
  template <typename T1> ndarray<T>& operator*=(const T1& x);
  template <typename T1> ndarray<T>& operator/=(const T1& x);

  template <typename T1> friend ndarray<T1> operator+(const ndarray<T1>& lhs, const ndarray<T1>& rhs);
  template <typename T1> friend ndarray<T1> operator-(const ndarray<T1>& lhs, const ndarray<T1>& rhs);

  // template <typename T1> friend ndarray<T> operator*(const ndarray<T>& lhs, const T1& rhs);
  template <typename T1> friend ndarray<T> operator*(const T1& lhs, const ndarray<T>& rhs) {return rhs * lhs;}
  template <typename T1> friend ndarray<T> operator/(const ndarray<T>& lhs, const T1& rhs);

  // element access
  T& operator[](size_t i) {return p[i];}
  const T& operator[](size_t i) const {return p[i];}

  template <typename F=float> // scalar multilinear interpolation
  T lerp(F x[]) const;

  template <int N, typename F=float> // vector multilinear interpolation
  T lerpv(F x[]) const;

  template <int N, typename F=float> // NxN tensor multilinear interpolation
  T lerpt(F x[]) const; 

  ndarray<T>& transpose(); // returns the ref to this
  ndarray<T> get_transpose() const; // only works for 2D arrays
  ndarray<T> get_transpose(const std::vector<size_t> order) const; // works for general tensors

  template <typename T1>
  void from_array(const ndarray<T1>& array1);

  void from_array(const T* p, const std::vector<size_t>& shape);

  void from_vector(const std::vector<T> &array);
  void copy_vector(const std::vector<T> &array);

  template <typename Iterator>
  void copy(Iterator first, Iterator last);

public: // subarray
  ndarray<T> subarray(const lattice&) const;

public: // file i/o; automatically determine format based on extensions
  static ndarray<T> from_file(const std::string& filename, const std::string varname="", diy::mpi::communicator comm = MPI_COMM_WORLD);
  bool read_file(const std::string& filename, const std::string varname="", diy::mpi::communicator comm = MPI_COMM_WORLD);
  bool to_file(const std::string& filename, const std::string varname="", diy::mpi::communicator comm = MPI_COMM_WORLD) const;

public: // i/o for binary file
  void read_binary_file(const std::string& filename) { ndarray_base::read_binary_file(filename); }
  void read_binary_file(FILE *fp);
  void read_binary_file_sequence(const std::string& pattern);
  void to_vector(std::vector<T> &out_vector) const;
  void to_binary_file(const std::string& filename) { ndarray_base::to_binary_file(filename); }
  void to_binary_file(FILE *fp);

  template <typename T1> void to_binary_file2(const std::string& f) const;

  void from_bov(const std::string& filename);
  void to_bov(const std::string& filename) const;

  void bil_add_block_raw(const std::string& filename, const std::vector<size_t>& SZ, const lattice& ext);

public: // i/o for vtk image data
  void to_vtk_image_data_file(const std::string& filename, const std::string varname=std::string()) const;
  void read_vtk_image_data_file_sequence(const std::string& pattern);
#if FTK_HAVE_VTK
  static int vtk_data_type();
  void from_vtu(vtkSmartPointer<vtkUnstructuredGrid> d, const std::string array_name=std::string());
  void from_vtk_image_data(vtkSmartPointer<vtkImageData> d, const std::string array_name=std::string()) { from_vtk_regular_data<>(d, array_name); }
  void from_vtk_array(vtkSmartPointer<vtkAbstractArray> d);
  void from_vtk_data_array(vtkSmartPointer<vtkDataArray> d);
  vtkSmartPointer<vtkImageData> to_vtk_image_data(std::string varname=std::string()) const; 
  vtkSmartPointer<vtkDataArray> to_vtk_data_array(std::string varname=std::string()) const; 
  
  template <typename VTK_REGULAR_DATA=vtkImageData> /*vtkImageData, vtkRectilinearGrid, or vtkStructuredGrid*/
  void from_vtk_regular_data(vtkSmartPointer<VTK_REGULAR_DATA> d, const std::string array_name=std::string());
#endif

public: // i/o for vtkStructuredGrid data
  void to_vtk_rectilinear_grid(const std::string& filename, const std::string varname=std::string()) const;

public: // i/o for hdf5
  static ndarray<T> from_h5(const std::string& filename, const std::string& name);
#if FTK_HAVE_HDF5
  bool read_h5_did(hid_t did);
  static hid_t h5_mem_type_id();
#endif

public: // i/o for parallel-netcdf
#if FTK_HAVE_PNETCDF
  void read_pnetcdf_all(int ncid, int varid, const MPI_Offset *st, const MPI_Offset *sz);
#endif

public: // i/o for adios2
  static ndarray<T> from_bp(const std::string& filename, const std::string& name, int step = -1, diy::mpi::communicator comm = MPI_COMM_WORLD);
  void read_bp(const std::string& filename, const std::string& varname, int step = -1, diy::mpi::communicator comm = MPI_COMM_WORLD) { ndarray_base::read_bp(filename, varname, step, comm); }

#if FTK_HAVE_ADIOS2
  void read_bp(
      adios2::IO &io, 
      adios2::Engine& reader,
      const std::string &varname, 
      int step = -1); // read all
  
  void read_bp(
      adios2::IO &io, 
      adios2::Engine& reader, 
      adios2::Variable<T>& var, 
      int step = -1);
#endif

public: // i/o for adios1
  static ndarray<T> from_bp_legacy(const std::string& filename, const std::string& varname, diy::mpi::communicator comm);
  bool read_bp_legacy(const std::string& filename, const std::string& varname, diy::mpi::communicator comm);
#if FTK_HAVE_ADIOS1
  bool read_bp_legacy(ADIOS_FILE *fp, const std::string& varname);
#endif

public: // i/o for png
  void read_png(const std::string& filename);
  void to_png(const std::string& filename) const;

public: // i/o for amira, see https://www.csc.kth.se/~weinkauf/notes/amiramesh.html
  bool read_amira(const std::string& filename); // T must be floatt

public: // pybind11
#if FTK_HAVE_PYBIND11
  ndarray(const pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast> &numpy_array);
  void from_numpy(const pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast> &numpy_array);
  pybind11::array_t<T, pybind11::array::c_style> to_numpy() const;
#endif
  void read_numpy(const std::string& filename);
  void to_numpy(const std::string& filename) const;

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

#if FTK_HAVE_MPI
  static MPI_Datatype mpi_datatype();
#endif
  static int nc_datatype();

  void copy_to_cuda_device();

public: // statistics & misc
  std::tuple<T, T> min_max() const;
  T maxabs() const;
  T resolution() const; // the min abs nonzero value
  
  ndarray<uint64_t> quantize() const; // quantization based on resolution

  ndarray<T> &perturb(T sigma); // add gaussian noise to the array
  ndarray<T> &clamp(T min, T max); // clamp data with min and max

private:
  std::vector<T> p;
    
#if 0 // FTK_HAVE_CUDA
  // arrays on GPU
  size_t *d_dims = NULL, *d_prod = NULL;
  T *d_p = NULL;
#endif
};

//////////////////////////////////

template <typename T> int ndarray<T>::type() const { return NDARRAY_TYPE_UNKNOWN; }
template <> inline int ndarray<double>::type() const { return NDARRAY_TYPE_DOUBLE; }
template <> inline int ndarray<float>::type() const { return NDARRAY_TYPE_FLOAT; }
template <> inline int ndarray<int>::type() const { return NDARRAY_TYPE_INT; }

template <typename T>
unsigned int ndarray<T>::hash() const
{
  unsigned int h0 = murmurhash2(p.data(), sizeof(T)*p.size(), 0);
  unsigned int h1 = murmurhash2(dims.data(), sizeof(size_t)*dims.size(), h0);
  return h1;
}

template <typename T>
template <typename T1>
ndarray<T>& ndarray<T>::operator*=(const T1& x)
{
  for (auto i = 0; i < p.size(); i ++)
    p[i] *= x;
  return *this;
}

template <typename T>
template <typename T1>
ndarray<T>& ndarray<T>::operator/=(const T1& x)
{
  for (auto i = 0; i < p.size(); i ++)
    p[i] /= x;
  return *this;
}

template <typename T>
ndarray<T>& ndarray<T>::operator+=(const ndarray<T>& x)
{
  if (empty()) *this = x;
  else {
    assert(this->shape() == x.shape());
    for (auto i = 0; i < p.size(); i ++)
      p[i] += x.p[i];
  }
  return *this;
}

template <typename T>
ndarray<T> operator+(const ndarray<T>& lhs, const ndarray<T>& rhs)
{
  ndarray<T> array;
  array.reshape(lhs);

  for (auto i = 0; i < array.nelem(); i ++)
    array[i] = lhs[i] + rhs[i];
  return array;
}

template <typename T, typename T1>
ndarray<T> operator*(const ndarray<T>& lhs, const T1& rhs) 
{
  ndarray<T> array;
  array.reshape(lhs);
  for (auto i = 0; i < array.nelem(); i ++)
    array[i] = lhs[i] * rhs;
  return array;
}

template <typename T, typename T1>
ndarray<T> operator/(const ndarray<T>& lhs, const T1& rhs) 
{
  ndarray<T> array;
  array.reshape(lhs);
  for (auto i = 0; i < array.nelem(); i ++)
    array[i] = lhs[i] / rhs;
  return array;
}

template <typename T>
void ndarray<T>::fill(T v)
{
  std::fill(p.begin(), p.end(), v);
}

template <typename T>
void ndarray<T>::fill(const std::vector<T>& values)
{
  p = values;
}

template <typename T>
void ndarray<T>::to_vector(std::vector<T> &out_vector) const{
  out_vector = p; 
}

template <typename T>
template <typename T1>
void ndarray<T>::from_array(const ndarray<T1>& array1)
{
  reshape(array1.shape());
  for (auto i = 0; i < p.size(); i ++)
    p[i] = static_cast<T>(array1[i]);
  ncd = array1.multicomponents();
  tv = array1.has_time();
}

template <typename T>
void ndarray<T>::from_array(const T *x, const std::vector<size_t>& shape)
{
  reshape(shape);
  memcpy(&p[0], x, nelem() * sizeof(T));
}

template <typename T>
void ndarray<T>::from_vector(const std::vector<T> &in_vector){
  for (int i=0;i<nelem();++i)
    if (i<in_vector.size())
      p[i] = in_vector[i];
    else break;
}

template <typename T>
void ndarray<T>::copy_vector(const std::vector<T> &array)
{
  p = array;
  reshape({p.size()});
}

template <typename T>
template <typename Iterator>
void ndarray<T>::copy(Iterator first, Iterator last)
{
  p.clear();
  std::copy(first, last, std::back_inserter(p));
  reshape({p.size()});
}

template <typename T>
void ndarray<T>::swap(ndarray& x)
{
  dims.swap(x.dims);
  s.swap(x.s);
  p.swap(x.p);
  std::swap(x.ncd, ncd);
  std::swap(x.tv, tv);
}

template <typename T>
ndarray<T> ndarray<T>::subarray(const lattice& l0) const
{
  lattice l(l0);
  if (l0.nd_cuttable() < nd()) {
    for (int i = 0; i < ncd; i ++) {
      l.starts_.insert(l.starts_.begin(), 0);
      l.sizes_.insert(l.starts_.begin(), this->shape(i));
    }
  }

  ndarray<T> arr(l.sizes());
  for (auto i = 0; i < arr.nelem(); i ++) {
    auto idx = l.from_integer(i);
    arr[i] = at(idx);
  }
  
  arr.ncd = ncd;
  return arr;
}

template <typename T>
void ndarray<T>::bil_add_block_raw(const std::string& filename, 
    const std::vector<size_t>& SZ, 
    const lattice& ext)
{
#if FTK_HAVE_MPI
  reshape(ext.sizes());
  std::vector<int> domain, st, sz;
 
  for (int i = 0; i < nd(); i ++) {
    domain.push_back(SZ[i]);
    st.push_back(ext.start(i));
    sz.push_back(ext.size(i));
  }

  BIL_Add_block_raw(nd(), domain.data(), st.data(), sz.data(), filename.c_str(), mpi_datatype(), (void**)&p[0]);
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_MPI);
#endif
}

template <typename T>
void ndarray<T>::read_binary_file(FILE *fp)
{
  auto s = fread(&p[0], sizeof(T), nelem(), fp);
  if (s != nelem())
    warn(FTK_ERR_FILE_CANNOT_READ_EXPECTED_BYTES);
}

template <typename T>
void ndarray<T>::to_binary_file(FILE *fp)
{
  fwrite(&p[0], sizeof(T), nelem(), fp);
}

template <typename T>
template <typename T1>
void ndarray<T>::to_binary_file2(const std::string& f) const
{
  ndarray<T1> array; 
  array.template from_array<T>(*this);
  array.to_binary_file(f);
}

template <typename T>
void ndarray<T>::read_binary_file_sequence(const std::string& pattern)
{
  const auto filenames = ftk::glob(pattern);
  if (filenames.size() == 0) return;

  std::vector<size_t> mydims = dims;
  mydims[nd() - 1] = filenames.size();
  reshape(mydims);
 
  size_t npt = std::accumulate(dims.begin(), dims.end()-1, 1, std::multiplies<size_t>());
  for (int i = 0; i < filenames.size(); i ++) {
    // fprintf(stderr, "loading %s\n", filenames[i].c_str());
    FILE *fp = fopen(filenames[i].c_str(), "rb");
    fread(&p[npt*i], sizeof(T), npt, fp);
    fclose(fp);
  }
}

#ifdef FTK_HAVE_VTK
template<> inline int ndarray<char>::vtk_data_type() {return VTK_CHAR;}
template<> inline int ndarray<unsigned char>::vtk_data_type() {return VTK_UNSIGNED_CHAR;}
template<> inline int ndarray<short>::vtk_data_type() {return VTK_SHORT;}
template<> inline int ndarray<unsigned short>::vtk_data_type() {return VTK_UNSIGNED_SHORT;}
template<> inline int ndarray<int>::vtk_data_type() {return VTK_INT;}
template<> inline int ndarray<unsigned int>::vtk_data_type() {return VTK_UNSIGNED_INT;}
template<> inline int ndarray<long>::vtk_data_type() {return VTK_LONG;}
template<> inline int ndarray<unsigned long>::vtk_data_type() {return VTK_UNSIGNED_LONG;}
template<> inline int ndarray<float>::vtk_data_type() {return VTK_FLOAT;}
template<> inline int ndarray<double>::vtk_data_type() {return VTK_DOUBLE;}

template <typename T>
inline void ndarray<T>::from_vtk_array(vtkSmartPointer<vtkAbstractArray> d)
{
  vtkSmartPointer<vtkDataArray> da = vtkDataArray::SafeDownCast(d);
  from_vtk_data_array(da);
}

template<typename T>
inline void ndarray<T>::from_vtk_data_array(
    vtkSmartPointer<vtkDataArray> da)
{
  const int nc = da->GetNumberOfComponents(), 
            ne = da->GetNumberOfTuples();
  if (nc > 1) {
    reshape(nc, ne);
    set_multicomponents(1);
  } else {
    reshape(ne);
    set_multicomponents(0);
  }

  for (auto i = 0; i < ne; i ++) {
    double *tuple = da->GetTuple(i);
    for (auto j = 0; j < nc; j ++)
      p[i*nc+j] = tuple[j];
  }
}

template <typename T>
inline void ndarray<T>::from_vtu(
    vtkSmartPointer<vtkUnstructuredGrid> d, 
    const std::string array_name)
{
  vtkSmartPointer<vtkDataArray> da = d->GetPointData()->GetArray(array_name.c_str());
  if (!da) da = d->GetPointData()->GetArray(0);
  from_vtk_data_array(da);
}

template <typename T>
template <typename VTK_REGULAR_DATA>
inline void ndarray<T>::from_vtk_regular_data(
    vtkSmartPointer<VTK_REGULAR_DATA> d, 
    const std::string array_name)
{
  vtkSmartPointer<vtkDataArray> da = d->GetPointData()->GetArray(array_name.c_str());
  if (!da) da = d->GetPointData()->GetArray(0);
  
  const int nd = d->GetDataDimension(), 
            nc = da->GetNumberOfComponents();

  if (nd == 2) {
    if (nc == 1) reshape(d->GetDimensions()[0], d->GetDimensions()[1]); // scalar field
    else {
      reshape(nc, d->GetDimensions()[0], d->GetDimensions()[1]); // vector field
      ncd = 1; // multicomponent array
    };
  } else if (nd == 3) {
    if (nc == 1) reshape(d->GetDimensions()[0], d->GetDimensions()[1], d->GetDimensions()[2]); // scalar field
    else {
      reshape(nc, d->GetDimensions()[0], d->GetDimensions()[1], d->GetDimensions()[2]);
      ncd = 1; // multicomponent array
    }
  } else {
    fprintf(stderr, "[FTK] fatal error: unsupported data dimension %d.\n", nd);
    assert(false);
  }

  for (auto i = 0; i < da->GetNumberOfTuples(); i ++) {
    double *tuple = da->GetTuple(i);
    for (auto j = 0; j < nc; j ++)
      p[i*nc+j] = tuple[j];
  }
}


template<typename T>
inline void ndarray<T>::to_vtk_image_data_file(const std::string& filename, const std::string varname) const
{
  // fprintf(stderr, "to_vtk_image_data_file, ncd=%zu\n", ncd);
  vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkXMLImageDataWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( to_vtk_image_data(varname) );
  writer->Write();
}

template <typename T>
inline vtkSmartPointer<vtkDataArray> ndarray<T>::to_vtk_data_array(const std::string varname) const
{
  vtkSmartPointer<vtkDataArray> d = vtkDataArray::CreateDataArray(vtk_data_type());
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
  memcpy(d->GetVoidPointer(0), p.data(), sizeof(T) * p.size()); // nelem());
  return d;
}

template<typename T>
inline vtkSmartPointer<vtkImageData> ndarray<T>::to_vtk_image_data(std::string varname) const
{
  vtkSmartPointer<vtkImageData> d = vtkImageData::New();
  // fprintf(stderr, "to_vtk_image_data, ncd=%zu\n", ncd);
  if (ncd) { // multicomponent
    if (nd() == 3) d->SetDimensions(shape(1), shape(2), 1);
    else d->SetDimensions(shape(1), shape(2), shape(3)); // nd == 4

    if (varname.empty()) varname = "vector";
  } else {
    if (nd() == 2) d->SetDimensions(shape(0), shape(1), 1);
    else d->SetDimensions(shape(0), shape(1), shape(2)); // nd == 3
    
    if (varname.empty()) varname = "scalar";
  }
  // d->GetPointData()->AddArray(to_vtk_data_array(multicomponent));
  d->GetPointData()->SetScalars(to_vtk_data_array(varname));

  return d;
}

template<typename T>
inline void ndarray<T>::read_vtk_image_data_file_sequence(const std::string& pattern)
{
  const auto filenames = glob(pattern);
  if (filenames.size() == 0) return;
  p.clear();

  ndarray<T> array;
  for (int t = 0; t < filenames.size(); t ++) {
    array.read_vtk_image_data_file(filenames[t]);
    p.insert(p.end(), array.p.begin(), array.p.end());
  }

  auto dims = array.dims;
  dims.push_back(filenames.size());
  reshape(dims);
}
#else
template<typename T>
inline void ndarray<T>::read_vtk_image_data_file_sequence(const std::string& pattern)
{
  fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
}

template<typename T>
inline void ndarray<T>::to_vtk_image_data_file(const std::string& filename, const std::string) const 
{
  fatal(FTK_ERR_NOT_BUILT_WITH_VTK);
}
#endif

template <typename T>
inline void ndarray<T>::read_netcdf(const std::string& filename, const std::string& varname, diy::mpi::communicator comm)
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

template <>
inline void ndarray<int>::read_netcdf(int ncid, int varid, int ndims, const size_t starts[], const size_t sizes[], diy::mpi::communicator comm)
{
#if FTK_HAVE_NETCDF
  std::vector<size_t> mysizes(sizes, sizes+ndims);
  std::reverse(mysizes.begin(), mysizes.end());
  reshape(mysizes);

  NC_SAFE_CALL( nc_get_vara_int(ncid, varid, starts, sizes, &p[0]) );
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

template <>
inline void ndarray<float>::read_netcdf(int ncid, int varid, int ndims, const size_t starts[], const size_t sizes[], diy::mpi::communicator comm)
{
#if FTK_HAVE_NETCDF
  std::vector<size_t> mysizes(sizes, sizes+ndims);
  std::reverse(mysizes.begin(), mysizes.end());
  reshape(mysizes);

#if NC_HAS_PARALLEL
  // NC_SAFE_CALL( nc_var_par_access(ncid, varid, NC_COLLECTIVE) );
  nc_var_par_access(ncid, varid, NC_COLLECTIVE);
#endif

  NC_SAFE_CALL( nc_get_vara_float(ncid, varid, starts, sizes, &p[0]) );
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

template <>
inline void ndarray<double>::read_netcdf(int ncid, int varid, int ndims, const size_t starts[], const size_t sizes[], diy::mpi::communicator comm)
{
#ifdef FTK_HAVE_NETCDF
  std::vector<size_t> mysizes(sizes, sizes+ndims);
  std::reverse(mysizes.begin(), mysizes.end());
  reshape(mysizes);

#if NC_HAS_PARALLEL
  // NC_SAFE_CALL( nc_var_par_access(ncid, varid, NC_COLLECTIVE) );
  nc_var_par_access(ncid, varid, NC_COLLECTIVE);
#endif

  NC_SAFE_CALL( nc_get_vara_double(ncid, varid, starts, sizes, &p[0]) );
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

template <>
inline void ndarray<double>::to_netcdf(int ncid, int varid, const size_t st[], const size_t sz[]) const
{
#ifdef FTK_HAVE_NETCDF
  fprintf(stderr, "st=%zu, %zu, %zu, %zu, sz=%zu, %zu, %zu, %zu\n", 
      st[0], st[1], st[2], st[3], sz[0], sz[1], sz[2], sz[3]);
  NC_SAFE_CALL( nc_put_vara_double(ncid, varid, st, sz, &p[0]) );
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

template <>
inline void ndarray<float>::to_netcdf(int ncid, int varid, const size_t starts[], const size_t sizes[]) const
{
#ifdef FTK_HAVE_NETCDF
  NC_SAFE_CALL( nc_put_vara_float(ncid, varid, starts, sizes, &p[0]) );
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

template <>
inline void ndarray<int>::to_netcdf(int ncid, int varid, const size_t starts[], const size_t sizes[]) const
{
#ifdef FTK_HAVE_NETCDF
  NC_SAFE_CALL( nc_put_vara_int(ncid, varid, starts, sizes, &p[0]) );
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

template <typename T>
inline void ndarray<T>::to_netcdf(int ncid, int varid) const
{
  std::vector<size_t> starts(dims.size(), 0), sizes(dims);
  std::reverse(sizes.begin(), sizes.end());

  to_netcdf(ncid, varid, &starts[0], &sizes[0]);
}

template <typename T>
inline void ndarray<T>::to_netcdf_multivariate(int ncid, int varids[]) const
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

template <typename T>
inline void ndarray<T>::to_netcdf_unlimited_time(int ncid, int varid) const
{
  std::vector<size_t> starts(dims.size()+1, 0), sizes(dims);
  sizes.push_back(1);
  std::reverse(sizes.begin(), sizes.end());
 
  // fprintf(stderr, "starts={%zu, %zu, %zu}, sizes={%zu, %zu, %zu}\n", 
  //     starts[0], starts[1], starts[2], sizes[0], sizes[1], sizes[2]);

  to_netcdf(ncid, varid, &starts[0], &sizes[0]);
}

template <typename T>
inline void ndarray<T>::to_netcdf_multivariate_unlimited_time(int ncid, int varids[]) const
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

template <typename T>
inline void ndarray<T>::read_netcdf(int ncid, int varid, const size_t starts[], const size_t sizes[], diy::mpi::communicator comm)
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

template <typename T>
inline void ndarray<T>::read_netcdf(int ncid, int varid, diy::mpi::communicator comm)
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

template <typename T>
inline void ndarray<T>::read_netcdf(int ncid, const std::string& varname, diy::mpi::communicator comm)
{
#ifdef FTK_HAVE_NETCDF
  int varid;
  NC_SAFE_CALL( nc_inq_varid(ncid, varname.c_str(), &varid) );
  read_netcdf(ncid, varid, comm);
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

template <typename T>
inline void ndarray<T>::read_netcdf(int ncid, const std::string& varname, const size_t starts[], const size_t sizes[], diy::mpi::communicator comm)
{
#ifdef FTK_HAVE_NETCDF
  int varid;
  NC_SAFE_CALL( nc_inq_varid(ncid, varname.c_str(), &varid) );
  read_netcdf(ncid, varid, starts, sizes, comm);
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_NETCDF);
#endif
}

template <typename T>
inline void ndarray<T>::read_netcdf(const std::string& filename, const std::string& varname, const size_t starts[], const size_t sizes[], diy::mpi::communicator comm)
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

template <typename T>
ndarray<T>::ndarray(const T *a, const std::vector<size_t> &dims_)
{
  from_array(a, dims_);
#if 0
  dims = dims_;
  s.resize(dims.size());

  for (size_t i = 0; i < nd(); i ++)
    if (i == 0) s[i] = 1;
    else s[i] = s[i-1]*dims[i-1];

  p.assign(a, a + s[nd()-1]);
#endif
}
  
template <typename T> 
template <typename I> 
void ndarray<T>::reshape(const int ndims, const I sz[])
{
  std::vector<size_t> sizes(ndims);
  for (int i = 0; i < ndims; i ++)
    sizes[i] = sz[i];
  reshape(sizes);
}

template <typename T>
template <typename T1>
void ndarray<T>::reshape(const ndarray<T1>& array)
{
  reshape(array.shape());
}

template <typename T>
void ndarray<T>::reshape(const std::vector<size_t> &dims_)
{
  dims = dims_;
  s.resize(dims.size());

  if (dims.size() == 0)
    fatal(FTK_ERR_NDARRAY_RESHAPE_EMPTY);

  for (size_t i = 0; i < nd(); i ++)
    if (i == 0) s[i] = 1;
    else s[i] = s[i-1]*dims[i-1];

  p.resize(s[nd()-1]*dims[nd()-1]);
}

template <typename T>
void ndarray<T>::reshape(const std::vector<size_t> &dims, T val)
{
  reshape(dims);
  std::fill(p.begin(), p.end(), val);
}

template <typename T>
std::tuple<T, T> ndarray<T>::min_max() const {
  T min = std::numeric_limits<T>::max(), 
    max = std::numeric_limits<T>::min();

  for (size_t i = 0; i < nelem(); i ++) {
    min = std::min(min, at(i));
    max = std::max(max, at(i));
  }

  return std::make_tuple(min, max);
}

template <typename T>
T ndarray<T>::maxabs() const 
{
  T r = 0;
  for (size_t i = 0; i < nelem(); i ++)
    r = std::max(r, std::abs(p[i]));

  return r;
}

template <typename T>
T ndarray<T>::resolution() const {
  T r = std::numeric_limits<T>::max();

  for (size_t i = 0; i < nelem(); i ++)
    if (p[i] != T(0))
      r = std::min(r, std::abs(p[i]));

  return r;
}

#if FTK_HAVE_MPI
template <> inline MPI_Datatype ndarray<double>::mpi_datatype() { return MPI_DOUBLE; }
template <> inline MPI_Datatype ndarray<float>::mpi_datatype() { return MPI_FLOAT; }
template <> inline MPI_Datatype ndarray<int>::mpi_datatype() { return MPI_INT; }
#endif

#if FTK_HAVE_NETCDF
template <> inline int ndarray<double>::nc_datatype() { return NC_DOUBLE; }
template <> inline int ndarray<float>::nc_datatype() { return NC_FLOAT; }
template <> inline int ndarray<int>::nc_datatype() { return NC_INT; }
#endif

#if FTK_HAVE_ADIOS2
template <typename T>
inline void ndarray<T>::read_bp(adios2::IO &io, adios2::Engine &reader, adios2::Variable<T>& var, int step)
{
  if (var) {
    // std::cerr << var << std::endl;
    // std::cerr << var.Shape() << std::endl;

    if (step == NDARRAY_ADIOS2_STEPS_UNSPECIFIED) {
      // nothing to do
    } else if (step == NDARRAY_ADIOS2_STEPS_ALL) {
      const size_t nsteps = var.Steps();
      var.SetStepSelection({0, nsteps-1});
    } else 
      var.SetStepSelection({step, 1});

    std::vector<size_t> shape(var.Shape());

    if (shape.size()) { // array type
      std::vector<size_t> zeros(shape.size(), 0);
      // fprintf(stderr, "%zu, %zu\n", shape.size(), zeros.size());
      reshape(shape);

      var.SetSelection({zeros, shape}); // read everything
      reader.Get<T>(var, p);

      std::reverse(shape.begin(), shape.end()); // we use a different dimension ordering than adios..
      reshape(shape);
    } else { // scalar type
      reshape(1);
      reader.Get<T>(var, p);
    }
  } else {
    throw FTK_ERR_ADIOS2_VARIABLE_NOT_FOUND;
    // fatal(FTK_ERR_ADIOS2_VARIABLE_NOT_FOUND);
  }
}

template <typename T>
inline void ndarray<T>::read_bp(adios2::IO &io, adios2::Engine &reader, const std::string &varname, int step)
{
  auto var = io.template InquireVariable<T>(varname);
  read_bp(io, reader, var, step);
}
#endif

#if FTK_HAVE_ADIOS1
template <typename T>
bool ndarray<T>::read_bp_legacy(ADIOS_FILE *fp, const std::string& varname)
{
  warn("reading bp file with legacy ADIOS1 API..");
  ADIOS_VARINFO *avi = adios_inq_var(fp, varname.c_str());
  if (avi == NULL)
    throw FTK_ERR_ADIOS2_VARIABLE_NOT_FOUND;

  adios_inq_var_stat(fp, avi, 0, 0);
  adios_inq_var_blockinfo(fp, avi);
  adios_inq_var_meshinfo(fp, avi);
    
  int nt = 1;
  uint64_t st[4] = {0, 0, 0, 0}, sz[4] = {0, 0, 0, 0};
  std::vector<size_t> mydims;
  
  for (int i = 0; i < avi->ndim; i++) {
    st[i] = 0;
    sz[i] = avi->dims[i];
    nt = nt * sz[i];
    mydims.push_back(sz[i]);
  }
  // fprintf(stderr, "%d, %d, %d, %d\n", sz[0], sz[1], sz[2], sz[3]);
  
  if (!mydims.empty()) {
    std::reverse(mydims.begin(), mydims.end());
    reshape(mydims);

    ADIOS_SELECTION *sel = adios_selection_boundingbox(avi->ndim, st, sz);
    assert(sel->type == ADIOS_SELECTION_BOUNDINGBOX);

    adios_schedule_read_byid(fp, sel, avi->varid, 0, 1, &p[0]);
    int retval = adios_perform_reads(fp, 1);
    
    adios_selection_delete(sel);
    return true; // avi->ndim;
  } else {
    if (avi->type == adios_integer) { // TODO: other data types
      reshape({1});
      p[0] = *((int*)avi->value);
      return true;
    }
    else return false;
  }
}
#endif

template <typename T>
ndarray<T> ndarray<T>::from_bp_legacy(const std::string& filename, const std::string& varname, diy::mpi::communicator comm)
{
  ndarray<T> arr;
  arr.read_bp_legacy(filename, varname, comm);
  return arr;
}

template <typename T>
bool ndarray<T>::read_bp_legacy(const std::string& filename, const std::string& varname, diy::mpi::communicator comm)
{
#if FTK_HAVE_ADIOS1
  adios_read_init_method( ADIOS_READ_METHOD_BP, comm, "" );
  ADIOS_FILE *fp = adios_read_open_file(filename.c_str(), ADIOS_READ_METHOD_BP, comm); 
  // adios_read_bp_reset_dimension_order(fp, 0);

  bool succ = read_bp_legacy(fp, varname);
  
  adios_read_finalize_method (ADIOS_READ_METHOD_BP);
  adios_read_close(fp);
  return succ;
#else
  fatal(FTK_ERR_NOT_BUILT_WITH_ADIOS1);
  return false;
#endif
}


template <typename T>
inline void ndarray<T>::copy_to_cuda_device()
{
#if 0 // FTK_HAVE_CUDA
  if (d_dims == NULL)
    cudaMalloc((void**)&d_dims, sizeof(size_t) * dims.size());
  cudaMemcpy(d_dims, dims.data(), sizeof(size_t) * dims.size(), cudaMemcpyHostToDevice);

  if (d_prod == NULL)
    cudaMalloc((void**)&d_prod, sizeof(size_t) * s.size());
  cudaMemcpy(d_prod, s.data(), sizeof(size_t) * s.size(), cudaMemcpyHostToDevice);

  if (d_p == NULL)
    cudaMalloc((void**)&d_p, sizeof(T) * nelem());
  cudaMemcpy(d_p, p.data(), sizeof(T) * p.size(), cudaMemcpyHostToDevice);
#else
  fprintf(stderr, "[FTK] fatal: FTK not compiled with CUDA.\n");
  assert(false);
#endif
}

template <typename T>
ndarray<T> ndarray<T>::from_file(const std::string& filename, const std::string varname, diy::mpi::communicator comm)
{
  ndarray<T> array;
  array.read_file(filename, varname, comm);
  return array;
}

template <typename T>
bool ndarray<T>::read_file(const std::string& filename, const std::string varname, diy::mpi::communicator comm)
{
  if (!file_exists(filename)) {
    warn(FTK_ERR_FILE_NOT_FOUND);
    return false;
  }

  auto ext = file_extension(filename);
  if (ext == FILE_EXT_BP) read_bp(filename, varname, -1, comm); // TODO: step
  else if (ext == FILE_EXT_NETCDF) read_netcdf(filename, varname, comm);
  else if (ext == FILE_EXT_VTI) read_vtk_image_data_file(filename, varname);
  else if (ext == FILE_EXT_HDF5) read_h5(filename, varname);
  else fatal(FTK_ERR_FILE_UNRECOGNIZED_EXTENSION);

  return true; // TODO: return read_* results
}

template <typename T>
ndarray<T> ndarray<T>::from_bp(const std::string& filename, const std::string& name, int step, diy::mpi::communicator comm)
{
  ndarray<T> array;
  array.read_bp(filename, name, step, comm);
    
  return array;
}

template <typename T>
ndarray<T> ndarray<T>::from_h5(const std::string& filename, const std::string& name)
{
  ndarray<T> array;
  array.read_h5(filename, name);
  return array;
}

#if FTK_HAVE_HDF5
template <typename T>
inline bool ndarray<T>::read_h5_did(hid_t did)
{
  auto sid = H5Dget_space(did); // space id
  auto type = H5Sget_simple_extent_type(sid);

  if (type == H5S_SIMPLE) {
    const int h5ndims = H5Sget_simple_extent_ndims(sid);
    hsize_t h5dims[h5ndims];
    H5Sget_simple_extent_dims(sid, h5dims, NULL);
    
    std::vector<size_t> dims(h5ndims);
    for (auto i = 0; i < h5ndims; i ++)
      dims[i] = h5dims[i];
    std::reverse(dims.begin(), dims.end()); // we use a different dimension ordering than hdf5..
    reshape(dims);
    
    H5Dread(did, h5_mem_type_id(), H5S_ALL, H5S_ALL, H5P_DEFAULT, p.data());
  } else if (type == H5S_SCALAR) {
    reshape(1);
    H5Dread(did, h5_mem_type_id(), H5S_ALL, H5S_ALL, H5P_DEFAULT, p.data());
  } else 
    fatal("unsupported h5 extent type");

  return true;
}

template <> inline hid_t ndarray<double>::h5_mem_type_id() { return H5T_NATIVE_DOUBLE; }
template <> inline hid_t ndarray<float>::h5_mem_type_id() { return H5T_NATIVE_FLOAT; }
template <> inline hid_t ndarray<int>::h5_mem_type_id() { return H5T_NATIVE_INT; }
template <> inline hid_t ndarray<unsigned long>::h5_mem_type_id() { return H5T_NATIVE_ULONG; }
template <> inline hid_t ndarray<unsigned int>::h5_mem_type_id() { return H5T_NATIVE_UINT; }
template <> inline hid_t ndarray<unsigned char>::h5_mem_type_id() { return H5T_NATIVE_UCHAR; }
#endif

template <typename T>
inline ndarray<T> ndarray<T>::slice(const lattice& l) const
{
  ndarray<T> array(l);
  for (auto i = 0; i < l.n(); i ++) {
    auto idx = l.from_integer(i);
    array[i] = at(idx);
  }
  return array;
}

template <typename T>
inline ndarray<T> ndarray<T>::slice(const std::vector<size_t>& st, const std::vector<size_t>& sz) const
{
  return slice(lattice(st, sz));
}

template <typename T>
inline ndarray<T> ndarray<T>::slice_time(size_t t) const 
{
  ndarray<T> array;
  std::vector<size_t> mydims(dims);
  mydims.resize(nd()-1);

  array.reshape(mydims);
  memcpy(&array[0], &p[t * s[nd()-1]], s[nd()-1] * sizeof(T));

  return array;
}

template <typename T>
inline std::vector<ndarray<T>> ndarray<T>::slice_time() const
{
  std::vector<ndarray<T>> arrays;
  const size_t nt = shape(nd()-1);
  for (size_t i = 0; i < nt; i ++) 
    arrays.push_back(slice_time(i));
  return arrays;
}

template <typename T>
std::ostream& ndarray<T>::print_shape(std::ostream& os) const
{
  os << "nd=" << nd() << ", array_dims={";
  for (size_t i = 0; i < dims.size(); i ++) 
    if (i < dims.size()-1) os << dims[i] << ", ";
    else os << dims[i] << "}, ";
  
  os << "prod={";
  for (size_t i = 0; i < s.size(); i ++) 
    if (i < s.size()-1) os << s[i] << ", ";
    else os << s[i] << "}, ";
  
  os << "size=" << p.size(); // << std::endl;
  return os;
}

template <typename T>
std::ostream& ndarray<T>::print(std::ostream& os) const
{
  print_shape(os);

  if (nd() == 1) {
    os << "[";
    for (size_t i = 0; i < dims[0]; i ++)
      if (i < dims[0]-1) os << at(i) << ", ";
      else os << at(i) << "]";
  } else if (nd() == 2) {
    os << "[";
    for (size_t j = 0; j < dims[1]; j ++) {
      os << "[";
      for (size_t i = 0; i < dims[0]; i ++)
        if (i < dims[0]-1) os << at(i, j) << ", ";
        else os << at(i, j) << "]";
      if (j < dims[1]-1) os << "], ";
      else os << "]";
    }
  } else if (nd() == 3) {
    os << "[";
    for (size_t k = 0; k < dims[2]; k ++) {
      for (size_t j = 0; j < dims[1]; j ++) {
        os << "[";
        for (size_t i = 0; i < dims[0]; i ++)
          if (i < dims[0]-1) os << at(i, j) << ", ";
          else os << at(i, j) << "]";
        if (j < dims[1]-1) os << "], ";
        else os << "]";
      }
      if (k < dims[2]-1) os << "], ";
      else os << "]";
    }
  }

  return os;
}

template <typename T>
ndarray<T>& ndarray<T>::transpose()
{
  *this = get_transpose();
  return *this;
}

template <typename T>
ndarray<T> ndarray<T>::get_transpose() const 
{
  ndarray<T> a;
  if (nd() == 0) return a;
  else if (nd() == 1) return *this;
  else if (nd() == 2) {
    a.reshape(dim(1), dim(0));
    for (auto i = 0; i < dim(0); i ++) 
      for (auto j = 0; j < dim(1); j ++)
        a(j, i) = at(i, j);
    return a;
  } else if (nd() == 3) {
    a.reshape(dim(2), dim(1), dim(0));
    for (auto i = 0; i < dim(0); i ++) 
      for (auto j = 0; j < dim(1); j ++)
        for (auto k = 0; k < dim(2); k ++)
          a(k, j, i) = at(i, j, k);
    return a;
  } else if (nd() == 4) {
    a.reshape(dim(3), dim(2), dim(1), dim(0));
    for (auto i = 0; i < dim(0); i ++) 
      for (auto j = 0; j < dim(1); j ++)
        for (auto k = 0; k < dim(2); k ++)
          for (auto l = 0; l < dim(3); l ++)
            a(l, k, j, i) = at(i, j, k, l);
    return a;
  } else {
    fatal(FTK_ERR_NDARRAY_UNSUPPORTED_DIMENSIONALITY);
    return a;
  }
}

template <typename T>
ndarray<T> ndarray<T>::concat(const std::vector<ndarray<T>>& arrays)
{
  ndarray<T> result;
  std::vector<size_t> result_shape = arrays[0].shape();
  result_shape.insert(result_shape.begin(), arrays.size());
  result.reshape(result_shape);

  const auto n = arrays[0].nelem();
  const auto n1 = arrays.size();

  for (auto i = 0; i < n; i ++)
    for (auto j = 0 ; j < n1; j ++)
      result[i*n1 + j] = arrays[j][i];

  return result;
}

template <typename T>
ndarray<T> ndarray<T>::stack(const std::vector<ndarray<T>>& arrays)
{
  ndarray<T> result;
  std::vector<size_t> result_shape = arrays[0].shape();
  result_shape.push_back(arrays.size());
  result.reshape(result_shape);

  const auto n = arrays[0].nelem();
  const auto n1 = arrays.size();

  for (auto j = 0 ; j < n1; j ++)
    for (auto i = 0; i < n; i ++)
      result[i + j*n] = arrays[j][i];

  return result;
}

#if FTK_HAVE_PYBIND11
template <typename T>
ndarray<T>::ndarray(const pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast> &numpy_array)
{
  from_numpy(numpy_array);
}
 
template <typename T>
void ndarray<T>::from_numpy(const pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast> &array)
{
  pybind11::buffer_info buf = array.request();
  std::vector<size_t> shape;
  for (auto i = 0; i < buf.ndim; i ++)
    shape.push_back(array.shape(i));
  reshape(shape);

  from_array((T*)buf.ptr, shape);
}

template <typename T>
pybind11::array_t<T, pybind11::array::c_style> ndarray<T>::to_numpy() const
{
  auto result = pybind11::array_t<T>(nelem());
  result.resize(shape());
  pybind11::buffer_info buf = result.request();

  T *ptr = (T*)buf.ptr;
  memcpy(ptr, data(), sizeof(T) * nelem());

  return result;
}
#endif

#if FTK_HAVE_PNG // TODO
template <typename T>
void ndarray<T>::to_png(const std::string& filename) const
{
  ndarray<unsigned char> buf;
  if (nd() == 2) 
    buf.reshape({4, dim(0), dim(1)});
  else if (nd() == 3 && dim(0) <= 4) 
    buf.reshape({4, dim(1), dim(2)});
  else 
    fatal("unable to save to png");

  // TODO
}
#else
template <typename T>
void ndarray<T>::read_png(const std::string& filename)
{
  fprintf(stderr, "[FTK] fatal error: FTK is not compiled with PNG.\n");
  assert(false);
}

template <typename T>
void ndarray<T>::to_png(const std::string& filename) const
{
  fprintf(stderr, "[FTK] fatal error: FTK is not compiled with PNG.\n");
  assert(false);
}
#endif

// see https://www.csc.kth.se/~weinkauf/notes/amiramesh.html
template <>
inline bool ndarray<float>::read_amira(const std::string& filename)
{
  auto find_and_jump = [](const char* buffer, const char* SearchString) {
    const char* FoundLoc = strstr(buffer, SearchString);
    if (FoundLoc) return FoundLoc + strlen(SearchString);
    return buffer;
  };

  FILE *fp = fopen(filename.c_str(), "rb"); 
  if (!fp) {
    warn(FTK_ERR_FILE_CANNOT_OPEN, filename);
    return false;
  }

  char buffer[2048];
  fread(buffer, sizeof(char), 2047, fp);
  buffer[2047] = '\0';

  if (!strstr(buffer, "# AmiraMesh BINARY-LITTLE-ENDIAN 2.1")) {
    warn(FTK_ERR_FILE_FORMAT_AMIRA, filename);
    fclose(fp);
    return false;
  }

  int xDim(0), yDim(0), zDim(0);
  sscanf(find_and_jump(buffer, "define Lattice"), "%d %d %d", &xDim, &yDim, &zDim);
  printf("\tAmriaMesh grid dimensions: %d %d %d\n", xDim, yDim, zDim);

  float xmin(1.0f), ymin(1.0f), zmin(1.0f);
  float xmax(-1.0f), ymax(-1.0f), zmax(-1.0f);

  sscanf(find_and_jump(buffer, "BoundingBox"), "%g %g %g %g %g %g", &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);
  printf("\tBoundingBox in x-Direction: [%g ... %g]\n", xmin, xmax);
  printf("\tBoundingBox in y-Direction: [%g ... %g]\n", ymin, ymax);
  printf("\tBoundingBox in z-Direction: [%g ... %g]\n", zmin, zmax);

  const bool bIsUniform = (strstr(buffer, "CoordType \"uniform\"") != NULL);
  printf("\tGridType: %s\n", bIsUniform ? "uniform" : "UNKNOWN");

  int NumComponents(0);
  if (strstr(buffer, "Lattice { float Data }"))
  { //Scalar field
    NumComponents = 1;
  } else {
    sscanf(find_and_jump(buffer, "Lattice { float["), "%d", &NumComponents);
  }
  printf("\tNumber of Components: %d\n", NumComponents);

  if (xDim <= 0 || yDim <= 0 || zDim <= 0
      || xmin > xmax || ymin > ymax || zmin > zmax
      || !bIsUniform || NumComponents <= 0)
  {
    warn(FTK_ERR_FILE_FORMAT_AMIRA);
    fclose(fp);
    return false;
  }

  const long idxStartData = strstr(buffer, "# Data section follows") - buffer;
  if (idxStartData > 0)
  {
    fseek(fp, idxStartData, SEEK_SET);
    fgets(buffer, 2047, fp);
    fgets(buffer, 2047, fp);

    reshape(NumComponents, xDim, yDim, zDim);
    set_multicomponents();
    read_binary_file(fp);
  }

  fclose(fp);
  return 0;
}

template <typename T>
ndarray<T>& ndarray<T>::perturb(T sigma)
{
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> d{0, sigma};

  for (auto i = 0; i < nelem(); i ++)
    p[i] = p[i] + d(gen);

  return *this;
}

template <typename T>
ndarray<T>& ndarray<T>::clamp(T min, T max)
{
  for (auto i = 0; i < nelem(); i ++)
    p[i] = ftk::clamp(p[i], min, max);

  return *this;
}

} // namespace ftk

///////// serialization

namespace diy {
  template <typename T> struct Serialization<ftk::ndarray<T>> {
    static void save(diy::BinaryBuffer& bb, const ftk::ndarray<T>& a) {
      diy::save(bb, a.dims);
      diy::save(bb, a.s);
      diy::save(bb, a.ncd);
      diy::save(bb, a.tv);
      diy::save(bb, a.p);
    }
   
    static void load(diy::BinaryBuffer& bb, ftk::ndarray<T>& a) {
      diy::load(bb, a.dims);
      diy::load(bb, a.s);
      diy::load(bb, a.ncd);
      diy::load(bb, a.tv);
      diy::load(bb, a.p);
    }
  };
} // namespace diy

#endif
