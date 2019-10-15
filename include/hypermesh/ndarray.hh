#ifndef _HYPERMESH_ARRAY_HH
#define _HYPERMESH_ARRAY_HH

#include <ftk/ftk_config.hh>
#include <vector>
#include <array>
#include <numeric>
#include <tuple>
#include <glob.h>

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

namespace hypermesh {

template <typename T>
struct ndarray {
  ndarray() {}
  ndarray(const std::vector<size_t> &dims) {reshape(dims);}
  ndarray(const T *a, const std::vector<size_t> &shape);

  size_t nd() const {return dims.size();}
  size_t dim(size_t i) const {return dims[i];}
  size_t shape(size_t i) const {return dim(i);}
  size_t nelem() const {return std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>());}
  std::vector<size_t> shape() const {return dims;}

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

private:
  std::vector<size_t> dims, s;
  std::vector<T> p;
};

template <typename T>
void ndarray<T>::to_vector(std::vector<T> &out_vector) const{
  out_vector = p; 
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
  fprintf(stderr, "[FTK] fatal error: your code is not compiled with netcdf.\n");
}

template <typename T>
inline void ndarray<T>::from_netcdf(int ncid, const std::string& varname, const size_t starts[], const size_t sizes[])
{
  fprintf(stderr, "[FTK] fatal error: your code is not compiled with netcdf.\n");
}

template <typename T>
inline void ndarray<T>::from_netcdf(const std::string& filename, const std::string& varname)
{
  fprintf(stderr, "[FTK] fatal error: your code is not compiled with netcdf.\n");
}

template <typename T>
inline void ndarray<T>::from_netcdf(const std::string& filename, const std::string& varname, const size_t starts[], const size_t sizes[])
{
  fprintf(stderr, "[FTK] fatal error: your code is not compiled with netcdf.\n");
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

}

#endif
