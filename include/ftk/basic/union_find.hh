#ifndef _UF_H
#define _UF_H

#include <vector>
#include <iostream>

// Implementation of weighted quick-union with path compression
// https://www.cs.princeton.edu/~rs/AlgsDS07/01UnionFind.pdf

namespace ftk {

template <class IdType=size_t>
struct quick_union
{
  quick_union(IdType size) {
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
struct weighted_quick_union : public quick_union<IdType>
{
  weighted_quick_union(IdType size) : quick_union<IdType>(size) {
    sz.resize(size, 1);
  }
  
  void unite(IdType p, IdType q) {
    IdType i = quick_union<IdType>::find(p);
    IdType j = quick_union<IdType>::find(q);

    if (sz[i] < sz[j]) {
      quick_union<IdType>::id[i] = j; 
      sz[j] += sz[i];
    } else {
      quick_union<IdType>::id[j] = i;
      sz[i] += sz[j];
    }
  }

private:
  std::vector<size_t> sz;
};

}

#endif
