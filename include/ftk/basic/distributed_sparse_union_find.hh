#ifndef _FTK_DISTRIBUTED_SPARSE_UNION_FIND_H
#define _FTK_DISTRIBUTED_SPARSE_UNION_FIND_H

#include <vector>
#include <iostream>
#include <map>

// Implementation of weighted quick-union with path compression
// https://www.cs.princeton.edu/~rs/AlgsDS07/01UnionFind.pdf

// Add the sparse representation by using Hash map/table 

namespace ftk {
template <class IdType=std::string>
struct distributed_sparse_union_find
{
  distributed_sparse_union_find() {
    
  }

  // Add and initialize elements
  void add(IdType i) {
    eles.insert(i); 
    id[i] = i;
    sz[i] = 1;
  }

  void set_parent(IdType i, IdType par) {
    id[i] = par; 
  }

  IdType parent(IdType i) {
    if(id.find(i) == id.end()) {
      return ""; 
    }

    return id[i]; 
  }

  bool is_root(IdType i) {
    if(!has(i)) {
      return false; 
    }

    return i == id[i]; 
  }

  bool has(IdType i) {
    return eles.find(i) != eles.end(); 
    // return id.find(i) != id.end(); 
  }

public:
  std::set<std::string> eles; 

private:
  // Use HashMap to support sparse union-find
  std::map<IdType, IdType> id;
  std::map<IdType, size_t> sz;
};

}

#endif
