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
#include <vtkDataReader.h>
#include <vtkNew.h>
#endif

namespace ftk {

template <typename T> struct ndarray;

// the non-template base class for ndarray
struct ndarray_base {
 
  virtual size_t size() const = 0;
  virtual bool empty() const = 0;

  size_t nd() const {return dims.size();}
  size_t dim(size_t i) const {return dims[i];}
  size_t shape(size_t i) const {return dim(i);}
  const std::vector<size_t> &shape() const {return dims;}
  size_t nelem() const {return std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>());}
  
  virtual void reshape(const std::vector<size_t> &dims_) = 0;
  void reshape(size_t ndims, const size_t sizes[]);
  void reshape(const ndarray_base& array); //! copy shape from another array
  template <typename T> void reshape(const ndarray<T>& array); //! copy shape from another array
  
  size_t index(const std::vector<size_t>& idx) const;
  size_t index(const std::vector<int>& idx) const;
  
  template <typename uint=size_t>
  std::vector<uint> from_index(uint i) const {return lattice().from_integer(i);}

  lattice get_lattice() const;

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

public: // vti i/o
  void read_vtk_image_data_file(const std::string& filename, const std::string array_name=std::string());
  virtual void read_vtk_image_data_file_sequence(const std::string& pattern) = 0;
#if FTK_HAVE_VTK
  virtual void from_vtk_image_data(vtkSmartPointer<vtkImageData> d, const std::string array_name=std::string()) = 0;
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

} // namespace ftk

#endif
