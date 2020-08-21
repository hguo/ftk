#ifndef _HYPERMESH_ARRAY_HH
#define _HYPERMESH_ARRAY_HH

#include <ftk/ftk_config.hh>
#include <ftk/object.hh>
#include <ftk/mesh/lattice.hh>
#include <vector>
#include <array>
#include <numeric>
#include <tuple>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <glob.h>

#if FTK_HAVE_ADIOS2
#include <adios2.h>
#endif

#if FTK_HAVE_CUDA
// #include <cuda.h>
#include <cuda_runtime.h>
#endif

#if FTK_HAVE_MPI
#include <mpi.h>
#endif

#if FTK_HAVE_HDF5
#include <hdf5.h>
#endif

#if FTK_HAVE_NETCDF
#include <netcdf.h>
#define NC_SAFE_CALL(call) {\
  int retval = call;\
  if (retval != 0) {\
    fprintf(stderr, "[NetCDF Error] %s, in file '%s', line %i.\n", nc_strerror(retval), __FILE__, __LINE__); \
    exit(EXIT_FAILURE); \
  }\
}
#endif

#if FTK_HAVE_PNETCDF
#include <pnetcdf.h>
#define PNC_SAFE_CALL(call) {\
  int retval = call;\
  if (retval != 0) {\
      fprintf(stderr, "[PNetCDF Error] %s, in file '%s', line %i.\n", ncmpi_strerror(retval), __FILE__, __LINE__); \
      exit(EXIT_FAILURE); \
  }\
}
#endif

#if FTK_HAVE_VTK
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkDataReader.h>
#include <vtkNew.h>
#endif

#if FTK_HAVE_PYBIND11
#include <pybind11/numpy.h>
#endif

namespace ftk {

template <typename T>
struct ndarray : object {
  ndarray() {}
  ndarray(const std::vector<size_t> &dims) {reshape(dims);}
  ndarray(const lattice& l) {reshape(l.sizes());}
  ndarray(const T *a, const std::vector<size_t> &shape);

  std::ostream& print(std::ostream& os) const;
  std::ostream& print_shape(std::ostream& os) const;

  size_t nd() const {return dims.size();}
  size_t dim(size_t i) const {return dims[i];}
  size_t shape(size_t i) const {return dim(i);}
  size_t nelem() const {return std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>());}
  bool empty() const  {return p.empty();}
  const std::vector<size_t> &shape() const {return dims;}

  lattice get_lattice() const;

  void fill(const std::vector<T>& values); //! fill values with std::vector
  void fill(const std::vector<std::vector<T>>& values); //! fill values

  const std::vector<T>& std_vector() const {return p;}

  const T* data() const {return p.data();}
  T* data() {return p.data();}

  void swap(ndarray& x);

  void reshape(const std::vector<size_t> &dims_);
  void reshape(const std::vector<size_t> &dims, T val);
  void reshape(size_t ndims, const size_t sizes[]);
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

  size_t index(const std::vector<size_t>& idx) const;
  size_t index(const std::vector<int>& idx) const;

  template <typename uint=size_t>
  std::vector<uint> from_index(uint i) const {return lattice().from_integer(i);}

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

  ndarray<T>& operator=(const ndarray<T>& x) {dims = x.dims; s = x.s; p = x.p; return *this;}
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

  ndarray<T> transpose() const; // only works for 2D arrays
  ndarray<T> transpose(const std::vector<size_t> order) const; // works for general tensors

  template <typename T1>
  void from_array(const ndarray<T1>& array1);

  void from_array(const T* p, const std::vector<size_t>& shape);

  void from_vector(const std::vector<T> &array);
  void copy_vector(const std::vector<T> &array);

  template <typename Iterator>
  void copy(Iterator first, Iterator last);

  void from_binary_file(const std::string& filename);
  void from_binary_file(FILE *fp);
  void from_binary_file_sequence(const std::string& pattern);
  void to_vector(std::vector<T> &out_vector) const;
  void to_binary_file(const std::string& filename);
  void to_binary_file(FILE *fp);

  template <typename T1> void to_binary_file2(const std::string& f) const;

  void from_bov(const std::string& filename);
  void to_bov(const std::string& filename) const;

  void from_vtk_image_data_file(const std::string& filename, const std::string array_name=std::string());
  void from_vtk_image_data_file_sequence(const std::string& pattern);
  void to_vtk_image_data_file(const std::string& filename, bool multicomponent=false) const;
#if FTK_HAVE_VTK
  static int vtk_data_type();
  void from_vtk_image_data(vtkSmartPointer<vtkImageData> d, const std::string array_name=std::string());
  vtkSmartPointer<vtkImageData> to_vtk_image_data(bool multicomponent=false) const;
  vtkSmartPointer<vtkDataArray> to_vtk_data_array(bool multicomponent=false) const;
#endif

  void from_h5(const std::string& filename, const std::string& name); 
#if FTK_HAVE_HDF5
  void from_h5(hid_t fid, const std::string& name);
  void from_h5(hid_t did);

  static hid_t h5_mem_type_id();
#endif

#if FTK_HAVE_PNETCDF
  void from_pnetcdf_all(int ncid, int varid, const MPI_Offset *st, const MPI_Offset *sz);
#endif

public: // adios2
#if FTK_HAVE_ADIOS2
  void from_adios2(
      adios2::IO &io, 
      adios2::Engine& reader, 
      const std::string &varname, 
      int step_start = 0);
#endif

public: // pybind11
#if FTK_HAVE_PYBIND11
  ndarray(const pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast> &numpy_array);
  void from_numpy(const pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast> &numpy_array);
  pybind11::array_t<T, pybind11::array::c_style> to_numpy() const;
#endif
  void from_numpy(const std::string& filename);
  void to_numpy(const std::string& filename) const;

public: // netcdf
  void from_netcdf(const std::string& filename, const std::string& varname, const size_t starts[], const size_t sizes[]);
  void from_netcdf(int ncid, const std::string& varname, const size_t starts[], const size_t sizes[]);
  void from_netcdf(int ncid, int varid, int ndims, const size_t starts[], const size_t sizes[]);
  void from_netcdf(int ncid, int varid, const size_t starts[], const size_t sizes[]);
  void from_netcdf(const std::string& filename, const std::string& varname);
  void from_netcdf(int ncid, const std::string& varname);
  void from_netcdf(int ncid, int varid);
  // void to_netcdf(int ncid, const std::string& varname);
  // void to_netcdf(int ncid, int varid);
  void to_netcdf(int ncid, int varid, const size_t starts[], const size_t sizes[]) const;
  void to_netcdf(int ncid, int varid) const;
  void to_netcdf_multivariate(int ncid, int varids[]) const;
  void to_netcdf_unlimited_time(int ncid, int varid) const;
  void to_netcdf_multivariate_unlimited_time(int ncid, int varids[]) const;

  static std::vector<std::string> glob(const std::string &pattern);

#if FTK_HAVE_MPI
  static MPI_Datatype mpi_datatype();
#endif
  static int nc_datatype();

  void copy_to_cuda_device();

  // statistics
  std::tuple<T, T> min_max() const;

private:
  std::vector<size_t> dims, s;
  std::vector<T> p;

#if FTK_HAVE_CUDA
  // arrays on GPU
  size_t *d_dims = NULL, *d_prod = NULL;
  T *d_p = NULL;
#endif
};

//////////////////////////////////

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

template <typename T, typename T1>
ndarray<T> operator*(const ndarray<T>& lhs, const T1& rhs) 
{
  ndarray<T> array;
  array.reshape(lhs);
  for (auto i = 0; i < array.nelem(); i ++)
    array[i] = lhs[i] * rhs;
  return array;
}

template <typename T>
lattice ndarray<T>::get_lattice() const {
  std::vector<size_t> st(nd(), 0), sz(dims);
  return lattice(st, sz);
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
}

template <typename T>
void ndarray<T>::from_binary_file(const std::string& filename)
{
  FILE *fp = fopen(filename.c_str(), "rb");
  from_binary_file(fp);
  fclose(fp);
}

template <typename T>
void ndarray<T>::from_binary_file(FILE *fp)
{
  fread(&p[0], sizeof(T), nelem(), fp);
}

template <typename T>
void ndarray<T>::to_binary_file(const std::string& f)
{
  FILE *fp = fopen(f.c_str(), "wb");
  to_binary_file(fp);
  fclose(fp);
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
std::vector<std::string> ndarray<T>::glob(const std::string& pattern)
{
  std::vector<std::string> filenames;
  glob_t results; 
  ::glob(pattern.c_str(), 0, NULL, &results); 
  for (int i=0; i<results.gl_pathc; i++)
    filenames.push_back(results.gl_pathv[i]); 
  globfree(&results);
  return filenames;
}

template <typename T>
void ndarray<T>::from_binary_file_sequence(const std::string& pattern)
{
  const auto filenames = ndarray<T>::glob(pattern);
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

template<typename T>
inline void ndarray<T>::to_vtk_image_data_file(const std::string& filename, bool multicomponent) const 
{
  vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkXMLImageDataWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData( to_vtk_image_data(multicomponent) );
  writer->Write();
}

template<typename T>
inline vtkSmartPointer<vtkDataArray> ndarray<T>::to_vtk_data_array(bool multicomponent) const
{
  vtkSmartPointer<vtkDataArray> d = vtkDataArray::CreateDataArray(vtk_data_type());
  if (multicomponent) d->SetName("vector");
  else d->SetName("scalar");
  if (multicomponent) {
    d->SetNumberOfComponents(shape(0));
    d->SetNumberOfTuples( std::accumulate(dims.begin()+1, dims.end(), 1, std::multiplies<size_t>()) );
  }
  else {
    d->SetNumberOfComponents(1);
    d->SetNumberOfTuples(nelem());
  }
  memcpy(d->GetVoidPointer(0), p.data(), sizeof(T) * p.size()); // nelem());
  return d;
}

template<typename T>
inline vtkSmartPointer<vtkImageData> ndarray<T>::to_vtk_image_data(bool multicomponent) const
{
  vtkSmartPointer<vtkImageData> d = vtkImageData::New();
  if (multicomponent) {
    if (nd() == 3) d->SetDimensions(shape(1), shape(2), 1);
    else d->SetDimensions(shape(1), shape(2), shape(3));
  } else {
    if (nd() == 2) d->SetDimensions(shape(0), shape(1), 1);
    else d->SetDimensions(shape(0), shape(1), shape(2));
  }
  // d->GetPointData()->AddArray(to_vtk_data_array(multicomponent));
  d->GetPointData()->SetScalars(to_vtk_data_array(multicomponent));

  return d;
}

template<typename T>
inline void ndarray<T>::from_vtk_image_data(
    vtkSmartPointer<vtkImageData> d, 
    const std::string array_name)
{
  vtkSmartPointer<vtkDataArray> da = d->GetPointData()->GetArray(array_name.c_str());
  if (!da) da = d->GetPointData()->GetArray(0);
  
  const int nd = d->GetDataDimension(), 
            nc = da->GetNumberOfComponents();

  if (nd == 2) {
    if (nc == 1) reshape(d->GetDimensions()[0], d->GetDimensions()[1]); // scalar field
    else reshape(nc, d->GetDimensions()[0], d->GetDimensions()[1]); // vector field
  } else if (nd == 3) {
    if (nc == 1) reshape(d->GetDimensions()[0], d->GetDimensions()[1], d->GetDimensions()[2]); // scalar field
    else reshape(nc, d->GetDimensions()[0], d->GetDimensions()[1], d->GetDimensions()[2]);
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
inline void ndarray<T>::from_vtk_image_data_file(const std::string& filename, const std::string array_name)
{
  vtkNew<vtkXMLImageDataReader> reader;
  reader->SetFileName(filename.c_str());
  reader->Update();
  from_vtk_image_data(reader->GetOutput(), array_name);
}

template<typename T>
inline void ndarray<T>::from_vtk_image_data_file_sequence(const std::string& pattern)
{
  const auto filenames = ndarray<T>::glob(pattern);
  if (filenames.size() == 0) return;
  p.clear();

  ndarray<T> array;
  for (int t = 0; t < filenames.size(); t ++) {
    array.from_vtk_image_data_file(filenames[t]);
    p.insert(p.end(), array.p.begin(), array.p.end());
  }

  auto dims = array.dims;
  dims.push_back(filenames.size());
  reshape(dims);
}
#else
template<typename T>
inline void ndarray<T>::from_vtk_image_data_file(const std::string& filename, const std::string array_name)
{
  fprintf(stderr, "[FTK] fatal error: FTK is not compiled with VTK.\n");
  assert(false);
}

template<typename T>
inline void ndarray<T>::from_vtk_image_data_file_sequence(const std::string& pattern)
{
  fprintf(stderr, "[FTK] fatal error: FTK is not compiled with VTK.\n");
  assert(false);
}

template<typename T>
inline void ndarray<T>::to_vtk_image_data_file(const std::string& filename, bool) const 
{
  fprintf(stderr, "[FTK] fatal error: FTK is not compiled with VTK.\n");
  assert(false);
}
#endif

#ifdef FTK_HAVE_NETCDF
template <>
inline void ndarray<float>::from_netcdf(int ncid, int varid, int ndims, const size_t starts[], const size_t sizes[])
{
  std::vector<size_t> mysizes(sizes, sizes+ndims);
  std::reverse(mysizes.begin(), mysizes.end());
  reshape(mysizes);

  NC_SAFE_CALL( nc_get_vara_float(ncid, varid, starts, sizes, &p[0]) );
}

template <>
inline void ndarray<double>::from_netcdf(int ncid, int varid, int ndims, const size_t starts[], const size_t sizes[])
{
  std::vector<size_t> mysizes(sizes, sizes+ndims);
  std::reverse(mysizes.begin(), mysizes.end());
  reshape(mysizes);

  NC_SAFE_CALL( nc_get_vara_double(ncid, varid, starts, sizes, &p[0]) );
}

template <>
inline void ndarray<double>::to_netcdf(int ncid, int varid, const size_t st[], const size_t sz[]) const
{
  fprintf(stderr, "st=%zu, %zu, %zu, %zu, sz=%zu, %zu, %zu, %zu\n", 
      st[0], st[1], st[2], st[3], sz[0], sz[1], sz[2], sz[3]);
  NC_SAFE_CALL( nc_put_vara_double(ncid, varid, st, sz, &p[0]) );
}

template <>
inline void ndarray<float>::to_netcdf(int ncid, int varid, const size_t starts[], const size_t sizes[]) const
{
  NC_SAFE_CALL( nc_put_vara_float(ncid, varid, starts, sizes, &p[0]) );
}

template <>
inline void ndarray<int>::to_netcdf(int ncid, int varid, const size_t starts[], const size_t sizes[]) const
{
  NC_SAFE_CALL( nc_put_vara_int(ncid, varid, starts, sizes, &p[0]) );
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
inline void ndarray<T>::from_netcdf(int ncid, int varid, const size_t starts[], const size_t sizes[])
{
  int ndims;
  NC_SAFE_CALL( nc_inq_varndims(ncid, varid, &ndims) );

  std::vector<size_t> mysizes(sizes, sizes+ndims);
  std::reverse(mysizes.begin(), mysizes.end());
  reshape(mysizes);

  from_netcdf(ncid, varid, ndims, starts, sizes);
}

template <typename T>
inline void ndarray<T>::from_netcdf(int ncid, int varid)
{
  int ndims;
  int dimids[4];
  size_t starts[4] = {0}, sizes[4] = {0};

  NC_SAFE_CALL( nc_inq_varndims(ncid, varid, &ndims) );
  NC_SAFE_CALL( nc_inq_vardimid(ncid, varid, dimids) );

  for (int i = 0; i < ndims; i ++)
    NC_SAFE_CALL( nc_inq_dimlen(ncid, dimids[i], &sizes[i]) );
  
  from_netcdf(ncid, varid, starts, sizes);
}

template <typename T>
inline void ndarray<T>::from_netcdf(int ncid, const std::string& varname)
{
  int varid;
  NC_SAFE_CALL( nc_inq_varid(ncid, varname.c_str(), &varid) );
  from_netcdf(ncid, varid);
}

template <typename T>
inline void ndarray<T>::from_netcdf(const std::string& filename, const std::string& varname)
{
  int ncid, varid;
  NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
  NC_SAFE_CALL( nc_inq_varid(ncid, varname.c_str(), &varid) );
  from_netcdf(ncid, varid);
  NC_SAFE_CALL( nc_close(ncid) );
}

template <typename T>
inline void ndarray<T>::from_netcdf(int ncid, const std::string& varname, const size_t starts[], const size_t sizes[])
{
  int varid;
  NC_SAFE_CALL( nc_inq_varid(ncid, varname.c_str(), &varid) );
  from_netcdf(ncid, varid, starts, sizes);
}

template <typename T>
inline void ndarray<T>::from_netcdf(const std::string& filename, const std::string& varname, const size_t starts[], const size_t sizes[])
{
  int ncid, varid;
  NC_SAFE_CALL( nc_open(filename.c_str(), NC_NOWRITE, &ncid) );
  NC_SAFE_CALL( nc_inq_varid(ncid, varname.c_str(), &varid) );
  from_netcdf(ncid, varid, starts, sizes);
  NC_SAFE_CALL( nc_close(ncid) );
}
#else
template <typename T>
inline void ndarray<T>::from_netcdf(int ncid, int varid, const size_t starts[], const size_t sizes[])
{
  fprintf(stderr, "[FTK] fatal error: FTK is not compiled with netcdf.\n");
  assert(false);
}

template <typename T>
inline void ndarray<T>::from_netcdf(int ncid, const std::string& varname, const size_t starts[], const size_t sizes[])
{
  fprintf(stderr, "[FTK] fatal error: FTK is not compiled with netcdf.\n");
  assert(false);
}

template <typename T>
inline void ndarray<T>::from_netcdf(const std::string& filename, const std::string& varname)
{
  fprintf(stderr, "[FTK] fatal error: FTK is not compiled with netcdf.\n");
  assert(false);
}

template <typename T>
inline void ndarray<T>::from_netcdf(const std::string& filename, const std::string& varname, const size_t starts[], const size_t sizes[])
{
  fprintf(stderr, "[FTK] fatal error: FTK is not compiled with netcdf.\n");
  assert(false);
}
#endif

template <typename T>
void ndarray<T>::reshape(size_t ndims, const size_t dims[])
{
  std::vector<size_t> mydims(dims, dims+ndims);
  reshape(mydims);
}

template <typename T>
ndarray<T>::ndarray(const T *a, const std::vector<size_t> &dims_)
{
  dims = dims_;
  s.resize(dims.size());

  for (size_t i = 0; i < nd(); i ++)
    if (i == 0) s[i] = 1;
    else s[i] = s[i-1]*dims[i-1];

  p.assign(a, a + s[nd()-1]);
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
size_t ndarray<T>::index(const std::vector<size_t>& idx) const {
  size_t i(idx[0]);
  for (size_t j = 1; j < nd(); j ++)
    i += idx[j] * s[j];
  return i;
}

template <typename T>
size_t ndarray<T>::index(const std::vector<int>& idx) const {
  std::vector<size_t> myidx(idx.begin(), idx.end());
  return index(myidx);
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
inline void ndarray<T>::from_adios2(adios2::IO &io, adios2::Engine &reader, const std::string &varname, int step_start)
{
  auto var = io.template InquireVariable<T>(varname);
  if (var) {
    var.SetStepSelection({step_start, 1});
    reshape(var.Shape());
    // var.SetSelection({{0, 0}, {var.Shape()[0], var.Shape()[1]}});
    reader.Get<T>(var, p);
  }
}
#endif

template <typename T>
inline void ndarray<T>::copy_to_cuda_device()
{
#if FTK_HAVE_CUDA
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

#if FTK_HAVE_HDF5
template <typename T>
inline void ndarray<T>::from_h5(const std::string& filename, const std::string& name)
{
  auto fid = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  from_h5(fid, name);
  H5Fclose(fid);
}

template <typename T>
inline void ndarray<T>::from_h5(hid_t fid, const std::string& name)
{
  auto did = H5Dopen2(fid, name.c_str(), H5P_DEFAULT);
  from_h5(did);
  H5Dclose(did);
}

template <typename T>
inline void ndarray<T>::from_h5(hid_t did)
{
  auto sid = H5Dget_space(did);
  const int h5ndims = H5Sget_simple_extent_ndims(sid);
  hsize_t h5dims[h5ndims];
  H5Sget_simple_extent_dims(sid, h5dims, NULL);
  
  std::vector<size_t> dims(h5ndims);
  for (auto i = 0; i < h5ndims; i ++)
    dims[i] = h5dims[i];
  std::reverse(dims.begin(), dims.end());
  reshape(dims);
  
  H5Dread(did, h5_mem_type_id(), H5S_ALL, H5S_ALL, H5P_DEFAULT, p.data());
}

template <> inline hid_t ndarray<double>::h5_mem_type_id() { return H5T_NATIVE_DOUBLE; }
template <> inline hid_t ndarray<float>::h5_mem_type_id() { return H5T_NATIVE_FLOAT; }
template <> inline hid_t ndarray<int>::h5_mem_type_id() { return H5T_NATIVE_INT; }
#else
template <typename T>
inline void ndarray<T>::from_h5(const std::string& filename, const std::string& name)
{
  fprintf(stderr, "[FTK] fatal: FTK not compiled with HDF5.\n");
  assert(false);
}
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
  os << "array_dims={";
  for (size_t i = 0; i < dims.size(); i ++) 
    if (i < dims.size()-1) os << dims[i] << ", ";
    else os << dims[i] << "}";
  os << std::endl;
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
  }

  return os;
}

template <typename T>
ndarray<T> ndarray<T>::transpose() const 
{
  ndarray<T> a;
  a.reshape(dim(1), dim(0));
  for (auto i = 0; i < dim(0); i ++) 
    for (auto j = 0; j < dim(1); j ++)
      a(j, i) = at(i, j);
  return a;
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

}

#endif
