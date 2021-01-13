#ifndef _FTK_DUF_HH
#define _FTK_DUF_HH

#include <ftk/ftk_config.hh>
#include <map>

#if FTK_HAVE_TBB
#include <tbb/concurrent_hash_map.h>
#endif

namespace ftk {

template <typename T>
struct duf {
  void unite(T, T);
  T find(T) const;
  
  size_t size() const { return parents.size(); }
  bool exists(T i) const; //  { return parents.find(i) != parents.end(); }

  std::function<int(T)> pid = [](int){return 0;};

  // void exchange();

private:
#if FTK_HAVE_TBB
  mutable tbb::concurrent_hash_map<T, T> parents;
#else
  mutable std::map<T, T> parents; // pointer to the parent in the local process
  std::mutex mutex;
#endif
};

//////
template <typename T>
bool duf<T>::exists(T i) const
{
#if FTK_HAVE_TBB
  typename tbb::concurrent_hash_map<T, T>::const_accessor it;
  return parents.find(it, i);
#else
  std::lock_guard<std::mutex> guard(mutex);
  return parents.find(i) != parents.end();
#endif
}

template <typename T>
void duf<T>::unite(T i, T j)
{
  i = find(i); // i <-- root(i)
  j = find(j); // j <-- root(j)

  if (i > j) std::swap(i, j); // ensure i<j

  // critical section
  {
#if !FTK_HAVE_TBB
    std::lock_guard<std::mutex> guard(mutex);
#endif
    parents.insert({j, i}); // parents[j] = i;
  }
}

template <typename T>
T duf<T>::find(T i) const
{
#if FTK_HAVE_TBB
  typename tbb::concurrent_hash_map<T, T>::accessor it;
  if (parents.find(it, i)) {
    while (i != it->second) {
      i = it->second;

      typename tbb::concurrent_hash_map<T, T>::const_accessor it1;
      parents.find(it1, i);
      it->second = it1->second;
    }
    return i;
  } else {
    parents.insert({i, i}); // parents[i] = i;
    return i;
  }
#else
  std::lock_guard<std::mutex> guard(mutex);
  
  if (parents.find(i) == parents.end()) { // if i does not exist, insert i to the graph and return i as the root
    parents.insert({i, i}); // parents[i] = i;
    return i;
  } else {
    while (i != parents[i]) {
      parents[i] = parents[parents[i]];
      i = parents[i];
    }
    return i;
  }
#endif
}

#if 0
template <typename T>
void duf<T>::exchange()
{
  std::set<T> local_roots;
  for (const auto &kv : parents) {
    if (kv.first == find(kv.second))
      local_roots.insert(kv.first);
  }
}
#endif

}

#endif
