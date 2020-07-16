#ifndef _FTK_NDARRAY_GROUP_HH
#define _FTK_NDARRAY_GROUP_HH

#include <ftk/ndarray.hh>

namespace ftk {

enum {
  NDARRAY_TYPE_FLOAT32,
  NDARRAY_TYPE_FLOAT64,
  NDARRAY_TYPE_INT32
  // TBA
};

struct data_group {
  static data_group* create() {return new data_group;}
  ~data_group() {
    for (auto kv : data) 
      this->free(kv.first);
  }

  std::map<std::string, std::pair<int/*type*/, void* /*ndarray*/>> data;

  template <typename T> const ndarray<T>& get(const std::string& key) const {
    static ndarray<T> local_null_array;
    const auto it = data.find(key);
    if (it == data.end()) {
      return local_null_array;
    } else {
      assert(type<T>() == it->second.first); // make sure type is consistent
      return *static_cast<ndarray<T>*>(it->second.second);
    }
  }

  template <typename T> void set(const std::string& key, const ndarray<T>& array) {
    auto it = data.find(key);
    if (it != data.end()) this->free(it->first);
    data[key] = std::make_pair(type<T>(), new ndarray<T>(array));
  }

  void free(const std::string& key) {
    auto it = data.find(key);
    if (it != data.end()) 
      free(it->second.first, it->second.second);
  }

  template <typename T> static int type();
  template <typename T> static ndarray<T>* cast(void* p) {return static_cast<ndarray<T>*>(p);}
  static void free(int type, void* p) {
    if (type == NDARRAY_TYPE_FLOAT32) delete cast<float>(p);
    else if (type == NDARRAY_TYPE_FLOAT64) delete cast<double>(p);
    else if (type == NDARRAY_TYPE_INT32) delete cast<int>(p);
    else assert(false);
  }

private: // non-copyable
  data_group() {};
  data_group(const data_group&) = delete;
  data_group& operator=(const data_group&) = delete;
};

/////
template<> int data_group::type<float>() {return NDARRAY_TYPE_FLOAT32;}
template<> int data_group::type<double>() {return NDARRAY_TYPE_FLOAT64;}
template<> int data_group::type<int>() {return NDARRAY_TYPE_INT32;}

}

#endif
