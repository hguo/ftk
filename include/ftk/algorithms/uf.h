#ifndef _UF_H
#define _UF_H

#include <vector>
#include <iostream>

// Implementation of weighted quick-union with path compression
// https://www.cs.princeton.edu/~rs/AlgsDS07/01UnionFind.pdf

template <class IdType, class NumberType>
struct UnionFind
{
  UnionFind(IdType size) {
    sz.resize(size, 1);
    id.resize(size);
    for (IdType i=0; i<size; i++)
      id[i] = i;
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
  IdType root(IdType i) {
    while (i != id[i]) {
      id[i] = id[id[i]];
      i = id[i];
    }
    return i;
  }
  
  std::vector<IdType> id;
  std::vector<NumberType> sz;
};

#endif
