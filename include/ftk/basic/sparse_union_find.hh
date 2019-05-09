#ifndef _FTK_SPARSE_UNION_FIND_H
#define _FTK_SPARSE_UNION_FIND_H

#include <vector>
#include <iostream>
#include <map>

// Implementation of weighted quick-union with path compression
// https://www.cs.princeton.edu/~rs/AlgsDS07/01UnionFind.pdf

namespace ftk {

template <class IdType=size_t>
struct sparse_quick_union
{
  sparse_quick_union(size_t _size) {
    size = _size; 
    reset();
  }

  void reset() {
    for (IdType i = 0; i < size; ++i)
      id[i] = i;
  }
  
  IdType root(IdType i) {
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
    id[i] = j;
  }

protected:
  // Use HashMap to support sparse union-find
  std::map<IdType, IdType> id;
  size_t size;
};

template <class IdType=size_t>
struct sparse_weighted_quick_union : public sparse_quick_union<IdType>
{
  sparse_weighted_quick_union(IdType size) : sparse_quick_union<IdType>(size) {
    sz.resize(size, 1);
  }
  
  void unite(IdType p, IdType q) {
    IdType i = sparse_quick_union<IdType>::root(p);
    IdType j = sparse_quick_union<IdType>::root(q);

    if (sz[i] < sz[j]) {
      sparse_quick_union<IdType>::id[i] = j; 
      sz[j] += sz[i];
    } else {
      sparse_quick_union<IdType>::id[j] = i;
      sz[i] += sz[j];
    }
  }

private:
  std::vector<size_t> sz;
};

}

#endif
