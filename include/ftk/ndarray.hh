#ifndef _HYPERMESH_ARRAY_HH
#define _HYPERMESH_ARRAY_HH

#include <ftk/ftk_config.hh>
#include <ftk/hypermesh/lattice.hh>
#include <vector>
#include <array>
#include <numeric>
#include <tuple>
#include <algorithm>
#include <cassert>
#include <glob.h>

#if FTK_HAVE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#if FTK_HAVE_MPI
#include <mpi.h>
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
#include <vtkImageData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkDataReader.h>
#endif

namespace ftk {

template <typename T>
struct ndarray {
  ndarray() {}
  ndarray(const std::vector<size_t> &dims) {reshape(dims);}
  ndarray(const T *a, const std::vector<size_t> &shape);

  size_t nd() const {return dims.size();}
  size_t dim(size_t i) const {return dims[i];}
  size_t shape(size_t i) const {return dim(i);}
  size_t nelem() const {return std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>());}
  bool empty() const  {return p.empty();}
  std::vector<size_t> shape() const {return dims;}

  lattice get_lattice() const;

  const T* data() const {return p.data();}
  T* data() {return p.data();}

  void swap(ndarray& x);

  void reshape(const std::vector<size_t> &dims_);
  void reshape(const std::vector<size_t> &dims, T val);
  void reshape(size_t ndims, const size_t sizes[]);

  void reshape(size_t n0) {reshape(std::vector<size_t>({n0}));}
  void reshape(size_t n0, size_t n1) {reshape({n0, n1});}
  void reshape(size_t n0, size_t n1, size_t n2) {reshape({n0, n1, n2});}
  void reshape(size_t n0, size_t n1, size_t n2, size_t n3) {reshape({n0, n1, n2, n3});}
  void reshape(size_t n0, size_t n1, size_t n2, size_t n3, size_t n4) {reshape({n0, n1, n2, n3, n4});}
  void reshape(size_t n0, size_t n1, size_t n2, size_t n3, size_t n4, size_t n5) {reshape({n0, n1, n2, n3, n4, n5});}
  void reshape(size_t n0, size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6) {reshape({n0, n1, n2, n3, n4, n5, n6});}

  size_t index(const std::vector<size_t>& idx) const;
  size_t index(const std::vector<int>& idx) const;

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

  T& operator[](size_t i) {return p[i];}
  const T& operator[](size_t i) const {return p[i];}

  template <typename T1>
  void from_array(const ndarray<T1>& array1);

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

  void from_vtk_image_data_file(const std::string& filename);
  void from_vtk_image_data_file_sequence(const std::string& pattern);
#if FTK_HAVE_VTK
  void from_vtk_image_data(vtkSmartPointer<vtkImageData> d);
#endif

#if FTK_HAVE_PNETCDF
  void from_pnetcdf_all(int ncid, int varid, const MPI_Offset *st, const MPI_Offset *sz);
#endif

  void from_netcdf(const std::string& filename, const std::string& varname, const size_t starts[], const size_t sizes[]);
  void from_netcdf(int ncid, const std::string& varname, const size_t starts[], const size_t sizes[]);
  void from_netcdf(int ncid, int varid, int ndims, const size_t starts[], const size_t sizes[]);
  void from_netcdf(int ncid, int varid, const size_t starts[], const size_t sizes[]);
  void from_netcdf(const std::string& filename, const std::string& varname);
  void from_netcdf(int ncid, const std::string& varname);
  void from_netcdf(int ncid, int varid);
  void to_netcdf(const std::string& filename, const std::string& varname);
  void to_netcdf(int ncid, const std::string& varname);
  void to_netcdf(int ncid, int varid);

  static std::vector<std::string> glob(const std::string &pattern);

  // statistics
  std::tuple<T, T> min_max() const;

#if FTK_HAVE_MPI
  static MPI_Datatype mpi_datatype();
#endif

#if FTK_HAVE_CUDA
  void copy_to_device();
#endif

private:
  std::vector<size_t> dims, s;
  std::vector<T> p;

#if FTK_HAVE_CUDA
  // arrays on GPU
  size_t *d_dims = NULL, *d_prod = NULL;
  T *d_p = NULL;
#endif
};

template <typename T>
lattice ndarray<T>::get_lattice() const {
  std::vector<size_t> st(nd(), 0), sz(dims);
  return lattice(st, sz);
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
template<typename T>
inline void ndarray<T>::from_vtk_image_data(vtkSmartPointer<vtkImageData> d)
{
  const int nd = d->GetDataDimension(), 
            nc = d->GetNumberOfScalarComponents();

  d->PrintSelf(std::cerr, vtkIndent(2));
  if (nd == 2) {
    if (nc == 1) { // scalar field
      reshape(d->GetDimensions()[0], d->GetDimensions()[1]);
      for (int j = 0; j < d->GetDimensions()[1]; j ++)
        for (int i = 0; i < d->GetDimensions()[0]; i ++)
          p.push_back(d->GetScalarComponentAsDouble(i, j, 0, 0));
    } else {
      reshape(d->GetNumberOfScalarComponents(), d->GetDimensions()[0], d->GetDimensions()[1]);
      for (int j = 0; j < d->GetDimensions()[1]; j ++)
        for (int i = 0; i < d->GetDimensions()[0]; i ++)
          for (int c = 0; c < d->GetNumberOfScalarComponents(); c ++)
            p.push_back(d->GetScalarComponentAsDouble(i, j, 0, c));
    }
  } else if (nd == 3) {
    if (nc == 1) { // scalar field
      reshape(d->GetDimensions()[0], d->GetDimensions()[1], d->GetDimensions()[2]);
      for (int k = 0; k < d->GetDimensions()[2]; k ++)
        for (int j = 0; j < d->GetDimensions()[1]; j ++)
          for (int i = 0; i < d->GetDimensions()[0]; i ++)
            p.push_back(d->GetScalarComponentAsDouble(i, j, 0, 0));
    } else {
      reshape(d->GetNumberOfScalarComponents(), d->GetDimensions()[0], d->GetDimensions()[1], d->GetDimensions()[2]);
      for (int k = 0; k < d->GetDimensions()[2]; k ++)
        for (int j = 0; j < d->GetDimensions()[1]; j ++)
          for (int i = 0; i < d->GetDimensions()[0]; i ++)
            for (int c = 0; c < d->GetNumberOfScalarComponents(); c ++)
              p.push_back(d->GetScalarComponentAsDouble(i, j, k, c));
    }
  } else {
    fprintf(stderr, "[FTK] fatal error: unsupported data dimension %d.\n", nd);
    assert(false);
  }
}

template<typename T>
inline void ndarray<T>::from_vtk_image_data_file(const std::string& filename)
{
  vtkNew<vtkXMLImageDataReader> reader;
  reader->SetFileName(filename.c_str());
  reader->Update();
  from_vtk_image_data(reader->GetOutput());
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
inline void ndarray<T>::from_vtk_image_data_file(const std::string& filename)
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
template <>
MPI_Datatype ndarray<double>::mpi_datatype()
{
  return MPI_DOUBLE;
}

template <>
MPI_Datatype ndarray<float>::mpi_datatype()
{
  return MPI_FLOAT;
}

template <>
MPI_Datatype ndarray<int>::mpi_datatype()
{
  return MPI_INT;
}
#endif

#if FTK_HAVE_CUDA
template <typename T>
void ndarray<T>::copy_to_device()
{
  if (d_dims == NULL)
    cudaMalloc((void**)&d_dims, sizeof(size_t) * dims.size());
  cudaMemcpy(d_dims, dims.data(), sizeof(size_t) * dims.size(), cudaMemcpyHostToDevice);

  if (d_prod == NULL)
    cudaMalloc((void**)&d_prod, sizeof(size_t) * s.size());
  cudaMemcpy(d_prod, s.data(), sizeof(size_t) * s.size(), cudaMemcpyHostToDevice);

  if (d_p == NULL)
    cudaMalloc((void**)&d_p, sizeof(T) * nelem());
  cudaMemcpy(d_p, p.data(), sizeof(T) * p.size(), cudaMemcpyHostToDevice);
}
#endif

}

#endif
