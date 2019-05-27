#ifndef _UF_H
#define _UF_H

#include <vector>
#include <iostream>

// Implementation of weighted quick-union with path compression
// https://www.cs.princeton.edu/~rs/AlgsDS07/01UnionFind.pdf

namespace ftk {

template <class IdType=size_t>
struct union_find
{
  union_find(IdType size) {
    id2parent.resize(size);
    reset();
  }

  // Initialization
  void reset() {
    for (IdType i=0; i<id2parent.size(); i++)
      id2parent[i] = i;
  }
  
  // Operations

  void unite(IdType p, IdType q) {
    IdType i = find(p);
    IdType j = find(q);
    id2parent[i] = j;
  }

  // Queries

  IdType find(IdType i) {
    while (i != id2parent[i]) {
      id2parent[i] = id2parent[id2parent[i]];
      i = id2parent[i];
    }
    return i;
  }
  
  bool same_set(IdType p, IdType q) {
    return find(p) == find(q);
  }

protected:
  std::vector<IdType> id2parent;
};

template <class IdType=size_t>
struct weighted_union_find : public union_find<IdType>
{
  weighted_union_find(IdType size) : union_find<IdType>(size) {
    sz.resize(size, 1);
  }
  
  // Operations

  void unite(IdType p, IdType q) {
    IdType i = union_find<IdType>::find(p);
    IdType j = union_find<IdType>::find(q);

    if (sz[i] < sz[j]) {
      union_find<IdType>::id2parent[i] = j; 
      sz[j] += sz[i];
    } else {
      union_find<IdType>::id2parent[j] = i;
      sz[i] += sz[j];
    }
  }

private:
  std::vector<size_t> sz;
};

}

#endif
