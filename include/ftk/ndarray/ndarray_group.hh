#ifndef _FTK_NDARRAY_GROUP_HH
#define _FTK_NDARRAY_GROUP_HH

#include <ftk/ndarray.hh>
#include <unordered_map>

namespace ftk {

struct ndarray_group : public std::unordered_map<std::string, ndarray_base*> {
  ndarray_group() {}
  virtual ~ndarray_group();

  void remove(const std::string key);

  bool has(const std::string key) const { return this->find(key) != this->end(); }

  void set(const std::string key, ndarray_base *ptr) { this->emplace(key, ptr); }
  template <typename T> void set(const std::string key, const ndarray<T> &arr);

  template <typename T> ndarray<T>* get_ptr(const std::string key) {
    // if (has(key)) return std::dynamic_pointer_cast<ndarray<T>>(at(key)); 
    if (has(key)) return (ndarray<T>*)at(key);
    else return nullptr;
  }

  template <typename T> ndarray<T> get(const std::string key) { return *get_ptr<T>(key); }

  // template <typename ... Args> ndarray_group(Args&&... args);
};

template <typename T>
void ndarray_group::set(const std::string key, const ndarray<T> &arr)
{
  // std::cerr << arr << std::endl;
  
  ndarray<T> *p = new ndarray<T>;
  ndarray_base *pb = p;
  *p = arr;

  // *(std::dynamic_pointer_cast<ndarray<T>>(parr)) = arr;
  // std::cerr << *(std::dynamic_pointer_cast<ndarray<T>>(parr)) << std::endl;
  this->set(key, pb);
}

inline void ndarray_group::remove(const std::string key)
{
  auto p = this->find(key);
  if (p != this->end()) {
    delete p->second;
    this->erase(p);
  }
}

inline ndarray_group::~ndarray_group()
{
  for (auto kv : *this) {
    delete kv.second;
    this->erase(kv.first);
  }
}

}

#endif
