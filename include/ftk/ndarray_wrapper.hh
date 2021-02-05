#ifndef _FTK_NDARRAY_WRAPPER_HH
#define _FTK_NDARRAY_WRAPPER_HH

#include <ftk/config.hh>
#include <ftk/ndarray.hh>
  
#define NDARRAY_TYPE_INVALID 0
#define NDARRAY_TYPE_FLOAT 1
#define NDARRAY_TYPE_DOUBLE 2
#define NDARRAY_TYPE_INT 3

#define FTK_NDARRAY_OP(rtn, call) {\
  if (type == NDARRAY_TYPE_FLOAT) rtn = ptr_float->call;\
  else if (type == NDARRAY_TYPE_DOUBLE) rtn = ptr_double->call;\
  else rtn = ptr_int->call;\
}

#define FTK_NDARRAY_RTN_OP(call) {\
  if (type == NDARRAY_TYPE_FLOAT) return ptr_float->call;\
  else if (type == NDARRAY_TYPE_DOUBLE) return ptr_double->call;\
  else return ptr_int->call;\
}

#define FTK_NDARRAY_VOID_OP(call) {\
  if (type == NDARRAY_TYPE_FLOAT) ptr_float->call;\
  else if (type == NDARRAY_TYPE_DOUBLE) ptr_double->call;\
  else ptr_int->call;\
}

namespace ftk {

// a non-template ndarray wrapper
struct ndarray_wrapper {
  explicit ndarray_wrapper(int t);
  ndarray_wrapper(const ndarray_wrapper& w); // shallow copy
  ndarray_wrapper& operator=(const ndarray_wrapper& w); // shallow copy
  ~ndarray_wrapper() {}

  void deep_copy(const ndarray_wrapper& w);

  template <typename T> std::shared_ptr<ndarray<T>> get() const;

public: // wrapping ndarray member funcs
  size_t nd() const { FTK_NDARRAY_RTN_OP(nd()); }
  size_t dim(size_t i) const { FTK_NDARRAY_RTN_OP( dim(i)); } 
  size_t size() const { FTK_NDARRAY_RTN_OP(size()); }
  bool empty() const { FTK_NDARRAY_RTN_OP(empty()); }
  const std::vector<size_t>& shape() const { FTK_NDARRAY_RTN_OP( shape() ); }

  void set_multicomponents(size_t c=1) { FTK_NDARRAY_VOID_OP( set_multicomponents(c) ); }
  bool multicomponents() const { FTK_NDARRAY_RTN_OP( multicomponents() ); }
  
  void set_has_time(bool b) { FTK_NDARRAY_VOID_OP( set_has_time(b) ); }
  bool has_time() const { FTK_NDARRAY_RTN_OP( has_time() ); }

  void reshape(const std::vector<size_t>& dims) { FTK_NDARRAY_VOID_OP( reshape(dims) ); }

private:
  int type = NDARRAY_TYPE_INVALID;
  std::shared_ptr<ndarray<float>> ptr_float;
  std::shared_ptr<ndarray<double>> ptr_double;
  std::shared_ptr<ndarray<int>> ptr_int;
};

//////////////
template <> std::shared_ptr<ndarray<float>> ndarray_wrapper::get<float>() const { return ptr_float; }
template <> std::shared_ptr<ndarray<double>> ndarray_wrapper::get<double>() const { return ptr_double; }
template <> std::shared_ptr<ndarray<int>> ndarray_wrapper::get<int>() const { return ptr_int; }

ndarray_wrapper::ndarray_wrapper(int t)
{
  type = t;
  if (t == NDARRAY_TYPE_FLOAT) ptr_float.reset(new ndarray<float>());
  else if (t == NDARRAY_TYPE_DOUBLE) ptr_double.reset(new ndarray<double>());
  else if (t == NDARRAY_TYPE_INT) ptr_int.reset(new ndarray<int>());
  else t = NDARRAY_TYPE_INVALID;
}

ndarray_wrapper::ndarray_wrapper(const ndarray_wrapper& w)
{
  *this = w;
}

ndarray_wrapper& ndarray_wrapper::operator=(const ndarray_wrapper& w)
{
  type = w.type;
  ptr_float = w.ptr_float;
  ptr_double = w.ptr_double;
  ptr_int = w.ptr_int;
  return *this;
}

}

#endif
