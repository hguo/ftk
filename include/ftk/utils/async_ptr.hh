#ifndef _FTK_ASYNC_PTR_HH
#define _FTK_ASYNC_PTR_HH

namespace ftk {

template <typename T>
struct async_ptr {
  template <typename Y> friend class async_ptr;

  async_ptr() : p(NULL), count(new unsigned int(0)) {}
  async_ptr(T *p_) : p(p_), count(new unsigned int(1)) {}
  template <typename Y> async_ptr(Y *p_) : p(p_), count(new unsigned int(1)) {}
  template <typename Y> async_ptr(async_ptr<Y> obj) : p(obj.p), count(obj.count) {
    if (obj.p != NULL)
      __sync_fetch_and_add(count, 1);
  }
  async_ptr(const async_ptr<T>& obj) : p(obj.p), count(obj.count) {
    if (obj.p != NULL) 
      __sync_fetch_and_add(count, 1);
  }
  async_ptr& operator=(const async_ptr<T>& obj) {
    p = obj.p;
    count = obj.count;
    if (obj.p != NULL)
      __sync_fetch_and_add(count, 1);
    return *this;
  }
#if 0
  async_ptr(async_ptr<T>&& obj) {
    p = obj.p;
    count = obj.count;
    obj.p = NULL;
    obj.count = NULL;
  }
  async_ptr& operator=(async_ptr<T>&& obj) {
    p = obj.p;
    count = obj.count;
    obj.p = NULL;
    obj.count = NULL;
    return *this;
  }
#endif
  ~async_ptr() { clear(); }

  void reset(T *p) {
    async_ptr<T>(p).swap(*this);
  }

  void swap(async_ptr<T>& obj) {
    std::swap(p, obj.p);
    std::swap(count, obj.count);
  }

  unsigned int get_count() const { return *count; }
  T* get() const { return *p; }
  T* operator->() const { return p; }
  T& operator*() const { return *p; }

private:
  void clear() {
    __sync_fetch_and_sub(count, 1);
    if (*count == 0) {
      if (p != NULL) delete p;
      delete count;
    }
  }

private:
  T *p = NULL;
  unsigned int *count;
};

}

#endif
