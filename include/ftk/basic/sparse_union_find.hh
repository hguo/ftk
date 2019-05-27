#ifndef _FTK_SPARSE_UNION_FIND_H
#define _FTK_SPARSE_UNION_FIND_H

#include <vector>
#include <map>
#include <set>
#include <iostream>

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
    // id2parent[i] = i; 
    id2parent.insert(std::make_pair (i, i)); 
    sz[i] = 1;
  }

  // Operations
  
  bool unite(IdType i, IdType j) {
    if(!has(i) || !has(j)) {
      return false ;
    }

    i = find(i);
    j = find(j);

    if (sz[i] < sz[j]) {
      // id2parent[i] = j; 
      id2parent.find(i)->second = j; 
      sz[j] += sz[i];
    } else {
      // id2parent[j] = i;
      id2parent.find(j)->second = i; 
      sz[i] += sz[j];
    }

    return true; 
  }

  // Queries

  bool has(IdType i) {
    return eles.find(i) != eles.end(); 
  }

  IdType parent(IdType i) {
    if(!has(i)) {
      return i; 
    }

    // return id2parent[i]; 
    return id2parent.find(i)->second; 
  }

  // Return the root of an element. 
  // Return empty string, if the element is not in the data structure 
  IdType find(IdType i) {
    if(!has(i)) {
      return i; 
    }

    IdType parent_i = parent(i); 
    while (i != parent_i) {
      // id2parent[i] = id2parent[id2parent[i]];
      id2parent.find(i)->second = parent(parent_i);

      i = parent_i;
      parent_i = parent(i); 
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

  std::vector<std::set<IdType>> get_sets() {
    std::map<IdType, std::set<IdType>> root2set; 
    for(auto ite = eles.begin(); ite != eles.end(); ++ite) {
      IdType root = find(*ite); 

      root2set[root].insert(*ite);
    }

    std::vector<std::set<IdType>> results; 
    for(auto ite = root2set.begin(); ite != root2set.end(); ++ite) {
      // if(is_root(*ite))
      results.push_back(ite->second); 
    }

    return results; 
  }

public:
  std::set<IdType> eles; 

private:
  // Use HashMap to support sparse union-find
  std::map<IdType, IdType> id2parent;
  std::map<IdType, size_t> sz;
};

}

#endif
