#ifndef _FTK_DUF_HH
#define _FTK_DUF_HH

#include <ftk/ftk_config.hh>
#include <map>

namespace ftk {

template <typename T>
struct duf {
  void unite(T, T);
  T find(T);

  std::function<int(T)> pid = [](int){return 0;};

  void exchange();

  std::map<T, T> parents; // pointer to the parent in the local process
  std::map<T, size_t> sz;
};

//////
template <typename T>
void duf<T>::unite(T i, T j)
{
  i = find(i); // i <-- root(i)
  j = find(j); // j <-- root(j)

  if (i > j) std::swap(i, j); // ensure i<j
  parents[j] = i;
}

template <typename T>
T duf<T>::find(T i)
{
  if (parents.find(i) == parents.end()) { // if i does not exist, insert i to the graph and return i as the root
    parents[i] = i;
    return i;
  } else {
    while (i != parents[i]) {
      parents[i] = parents[parents[i]];
      i = parents[i];
    }
    return i;
  }
}

template <typename T>
void duf<T>::exchange()
{
  std::set<T> local_roots;
  for (const auto &kv : parents) {
    if (kv.first == find(kv.second))
      local_roots.insert(kv.first);
  }
  
}

}

#endif
