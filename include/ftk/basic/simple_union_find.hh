#ifndef _FTK_SIMPLE_UNION_FIND_H
#define _FTK_SIMPLE_UNION_FIND_H

#include <vector>
#include <iostream>

// Reference 
  // Online slides: https://www.cs.princeton.edu/~rs/AlgsDS07/01UnionFind.pdf
  // Paper: "Worst-case Analysis of Set Union Algorithms"
  
namespace ftk {

template <class IdType=size_t>
struct simple_union_find
{
  simple_union_find(IdType size) : id2parent(), sz() {
    id2parent.resize(size);
    sz.resize(size, 1);

    reset();
  }

  // Initialization
  void reset() {
    for (IdType i=0; i<id2parent.size(); i++)
      id2parent[i] = i;
  }
  
  // Operations

  // Union by size
  bool unite(IdType p, IdType q) {
    IdType i = find(p);
    IdType j = find(q);

    if (sz[i] < sz[j]) {
      id2parent[i] = j; 
      sz[j] += sz[i];
    } else {
      id2parent[j] = i;
      sz[i] += sz[j];
    }

    return true; 
  }

  // Queries

  // Find the root of an element. 
    // Path compression by path halving method
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

private:
  std::vector<IdType> id2parent;
  std::vector<size_t> sz;
};

}

#endif
