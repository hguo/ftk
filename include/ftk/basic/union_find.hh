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
    id.resize(size);
    reset();
  }

  void reset() {
    for (IdType i=0; i<id.size(); i++)
      id[i] = i;
  }
  
  IdType find(IdType i) {
    while (i != id[i]) {
      id[i] = id[id[i]];
      i = id[i];
    }
    return i;
  }
  
  bool same_set(IdType p, IdType q) {
    return find(p) == find(q);
  }

  void unite(IdType p, IdType q) {
    IdType i = find(p);
    IdType j = find(q);
    id[i] = j;
  }

protected:
  std::vector<IdType> id;
};

template <class IdType=size_t>
struct weighted_union_find : public union_find<IdType>
{
  weighted_union_find(IdType size) : union_find<IdType>(size) {
    sz.resize(size, 1);
  }
  
  void unite(IdType p, IdType q) {
    IdType i = union_find<IdType>::find(p);
    IdType j = union_find<IdType>::find(q);

    if (sz[i] < sz[j]) {
      union_find<IdType>::id[i] = j; 
      sz[j] += sz[i];
    } else {
      union_find<IdType>::id[j] = i;
      sz[i] += sz[j];
    }
  }

private:
  std::vector<size_t> sz;
};

}

#endif
