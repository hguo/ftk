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

  // Initialization
  // Add and initialize elements
  void add(IdType i) {
    eles.insert(i); 
    id2parent[i] = i;
    sz[i] = 1;
  }

  // Operations
  
  void unite(IdType i, IdType j) {
    if(!has(i) || !has(j)) {
      return ;
    }

    i = find(i);
    j = find(j);

    if (sz[i] < sz[j]) {
      id2parent[i] = j; 
      sz[j] += sz[i];
    } else {
      id2parent[j] = i;
      sz[i] += sz[j];
    }
  }

  // Queries

  bool has(IdType i) {
    return eles.find(i) != eles.end(); 
  }

  IdType parent(IdType i) {
    if(!has(i)) {
      return ""; 
    }

    return id2parent[i]; 
  }

  // Return the root of an element. 
  // Return empty string, if the element is not in the data structure 
  IdType find(IdType i) {
    if(!has(i)) {
      return ""; 
    }

    while (i != id2parent[i]) {
      id2parent[i] = id2parent[id2parent[i]];
      i = id2parent[i];
    }

    return i; 
  }

  bool is_root(IdType i) {
    if(!has(i)) {
      return false; 
    }

    return i == id2parent[i]; 
  }

  bool same_set(IdType i, IdType j) {
    return find(i) == find(j);
  }

public:
  std::set<std::string> eles; 

private:
  // Use HashMap to support sparse union-find
  std::map<IdType, IdType> id2parent;
  std::map<IdType, size_t> sz;
};

}

#endif
