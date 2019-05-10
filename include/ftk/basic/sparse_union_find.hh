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
struct sparse_quick_union
{
  sparse_quick_union() {
    
  }

  IdType root(IdType i) {
    if(id.find(i) == id.end()) {
      id[i] = i;
      sz[i] = 1;

      return i; 
    }

    while (i != id[i]) {
      id[i] = id[id[i]];
      i = id[i];
    }

    return i; 
  }

  bool find(IdType p, IdType q) {
    return root(p) == root(q);
  }
  
  void unite(IdType p, IdType q) {
    IdType i = root(p);
    IdType j = root(q);

    if (sz[i] < sz[j]) {
      id[i] = j; 
      sz[j] += sz[i];
    } else {
      id[j] = i;
      sz[i] += sz[j];
    }
  }

private:
  // Use HashMap to support sparse union-find
  std::map<IdType, IdType> id;
  std::map<IdType, size_t> sz;
};

}

#endif
