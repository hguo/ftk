#ifndef _FTK_DUF_HH
#define _FTK_DUF_HH

#include <ftk/config.hh>
#include <ftk/utils/gather.hh>
#include <map>
#include <mutex>

namespace ftk {

template <typename T, class Compare = std::less<T>>
struct duf {
  duf(diy::mpi::communicator c = MPI_COMM_WORLD) : comm(c) {}

  void unite(T, T);
  T find(T) const;
  
  size_t size() const { return parents.size(); }
  bool has(T i) const; //  { return parents.find(i) != parents.end(); }

  void sync();
  // void print() const;

  std::set<T, Compare> get_roots() const;

private:
  mutable std::map<T, T, Compare> parents; // parents in the local process
  mutable std::mutex mutex;

  diy::mpi::communicator comm;
};

//////
template <typename T, class Compare>
bool duf<T, Compare>::has(T i) const
{
  std::lock_guard<std::mutex> guard(mutex);
  return parents.find(i) != parents.end();
}

template <typename T, class Compare>
void duf<T, Compare>::unite(T i, T j)
{
  i = find(i); // i <-- root(i)
  j = find(j); // j <-- root(j)
  
  if (i == j) return;
  if (i > j) std::swap(i, j); // ensure i<j

  // critical section
  {
    std::lock_guard<std::mutex> guard(mutex);
    parents[j] = i; 
  }
}

template <typename T, class Compare>
T duf<T, Compare>::find(T i) const
{
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
}

#if 0
template <typename T>
void duf<T>::print() const
{
  for (const auto &kv : parents) 
    fprintf(stderr, "rank=%d, val=%d, root=%d\n", 
        comm.rank(), kv.first, find(kv.second));
}
#endif

template <typename T, class Compare>
std::set<T, Compare> duf<T, Compare>::get_roots() const
{
  std::set<T, Compare> roots;
  for (const auto &kv : parents) {
    if (kv.first == find(kv.second))
      roots.insert(kv.first);
  }
  return roots;
}

template <typename T, class Compare>
void duf<T, Compare>::sync()
{
  if (comm.size() == 1) return; // no need to sync

#if 0 // WIP
  for (auto &kv : parents)
    kv.second = find(kv.second); // set to root

  fprintf(stderr, "before: rank=%d, #parents=%zu\n", comm.rank(), parents.size());

  diy::mpi::allgather<std::map<T, T, Compare>>(comm, parents, parents, 
      [](const std::map<T, T, Compare>& in, std::map<T, T, Compare>& out) {
        for (const auto &kv : in) {
          auto it = out.find( kv.first );
          if (it == out.end()) 
            out.insert( kv );
          else if (kv.second < it->second) // (Compare(kv.second, it->second))
          {
            // fprintf(stderr, "change..\n");
            it->second = kv.second;
          }
        }
      });
  
  fprintf(stderr, "after: rank=%d, #parents=%zu\n", comm.rank(), parents.size());
#endif

#if 0
  while (1) {
    std::set<T, Compare> local_roots = get_roots();
    std::set<T, Compare> remote_roots;
    diy::mpi::allgather(comm, local_roots, remote_roots);
    // fprintf(stderr, "rank=%d, #local_roots=%d, #remote_roots=%d\n", comm.rank(), local_roots.size(), remote_roots.size());

    std::map<T, T, Compare> local_updated_parents, remote_updated_parents;
    for (const auto v : remote_roots) {
      // if (has(v)) std::cerr << "rank=" << comm.rank() << " has remote root " << v << std::endl;
      // else std::cerr << "rank=" << comm.rank() << " does not have remote root " << v << std::endl;

      if (find(v) != v) {
        local_updated_parents[v] = find(v);
        // fprintf(stderr, "rank=%d, updating %d\n", comm.rank(), v);
        // std::cerr << "rank=" << comm.rank() << ", updating " << v << std::endl;
      } // else std::cerr << "remote root " << v << " cannot find ancestor" << std::endl;
    }
    
    diy::mpi::allgather(comm, local_updated_parents, remote_updated_parents);
    fprintf(stderr, "rank=%d, #updated_remote_parents=%d\n", comm.rank(), remote_updated_parents.size());
    for (const auto &kv : remote_updated_parents) {
      // unite(kv.first, kv.second);
      parents[kv.first] = kv.second;
    }

    if (remote_updated_parents.size() == 0)
      break;
  }
#endif
}

} // namespace ftk

#endif
