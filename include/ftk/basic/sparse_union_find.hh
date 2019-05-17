#ifndef _FTK_SPARSE_UNION_FIND_H
#define _FTK_SPARSE_UNION_FIND_H

#include <vector>
#include <iostream>
#include <map>

// Implementation of weighted quick-union with path compression
// https://www.cs.princeton.edu/~rs/AlgsDS07/01UnionFind.pdf

// Add the sparse representation by using Hash map/table 

namespace ftk {
template <class IdType=std::string>
struct sparse_union_find
{
  sparse_union_find() {
    
  }

  // Add and initialize elements
  void add(IdType i) {
    eles.push_back(i); 
    id[i] = i;
    sz[i] = 1;
  }

  void set_parent(IdType i, IdType par) {
    id[i] = par; 
  }

  // Return the root of an element. 
  // Return empty string, if the element is not in the data structure 
  IdType find(IdType i) {
    if(id.find(i) == id.end()) {
      return ""; 
    }

    while (i != id[i]) {
      id[i] = id[id[i]];
      i = id[i];
    }

    return i; 
  }

  IdType parent(IdType i) {
    if(id.find(i) == id.end()) {
      return ""; 
    }

    return id[i]; 
  }
  
  void unite(IdType i, IdType j) {
    i = find(i);
    j = find(j);

    if (sz[i] < sz[j]) {
      id[i] = j; 
      sz[j] += sz[i];
    } else {
      id[j] = i;
      sz[i] += sz[j];
    }
  }

  bool same_set(IdType i, IdType j) {
    return find(i) == find(j);
  }

  bool is_root(IdType i) {
    return i == id[i]; 
  }

  bool has(IdType i) {
    return id.find(i) != id.end(); 
  }

public:
  std::vector<std::string> eles; 

private:
  // Use HashMap to support sparse union-find
  std::map<IdType, IdType> id;
  std::map<IdType, size_t> sz;
};

}

#endif
