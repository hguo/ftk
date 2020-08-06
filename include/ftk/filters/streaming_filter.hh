#ifndef _FTK_STREAMING_FILTER_HH
#define _FTK_STREAMING_FILTER_HH

#include <deque>
#include <vector>
#include <functional>
#include <ftk/ndarray/conv.hh>

namespace ftk {

template <typename T/*data type such as ndarray*/, typename KT/*kernel type such as float/double*/>
struct streaming_filter {
  streaming_filter() : cursor(0) {}
  ~streaming_filter() {};

  void set_callback(std::function<void(int, const T&)> f) {callback = f;}

  void set_gaussian_kernel(KT sigma, int size) {set_kernel(gaussian_kernel(sigma, size));}
  void set_kernel(const std::vector<KT>&);

  const std::vector<KT>& get_kernel() const {return kernel;}
  int kernel_size() const {return kernel.size();}
  int half_kernel_size() const {return (kernel_size() + 1) / 2;}

  void push(const T& d);
  void finish();

protected:
  bool ready() const {return data.size() >= half_kernel_size();}
  T get() const;

private:
  std::vector<KT> kernel;
  std::deque<T> data;
  mutable int cursor;
  bool finishing = false;

  std::function<void(int, const T&)> callback;
  mutable int n_fetched_timesteps = 0;
};


/////////
template <typename T, typename KT>
void streaming_filter<T, KT>::set_kernel(const std::vector<KT>& k)
{
  kernel = k;
  assert(kernel.size() % 2 == 1);
}

template <typename T, typename KT>
void streaming_filter<T, KT>::push(const T& d)
{
  assert(kernel.size() % 2 == 1);

  data.push_back(d);
  if (data.size() > kernel.size()) {
    data.pop_front();
    cursor --;
  }

  if (callback) 
    if (ready()) {
      callback(n_fetched_timesteps, get());
    }
}

template <typename T, typename KT>
T streaming_filter<T, KT>::get() const
{
  assert(ready());

  // fprintf(stderr, "cursor=%d, #deque=%zu\n", cursor, data.size());
  T result;
  for (int i = 0; i < kernel_size(); i++) {
    int j;
    if (finishing) 
      j = std::min((int)data.size() - 1, i);
    else 
      j = std::max(0, i + cursor - half_kernel_size() + 1);
    // fprintf(stderr, "--j=%d\n", j);
    result += kernel[i] *  data[j];
  }

  if (finishing) cursor --;
  else cursor ++;

  n_fetched_timesteps ++;
  return result;
}

template <typename T, typename KT>
void streaming_filter<T, KT>::finish()
{
  finishing = true;
  while (1) {
    data.pop_front();
    if (ready() && callback) {
      callback(n_fetched_timesteps, get());
    }
    if (!ready()) break;
  }
}

}

#endif
