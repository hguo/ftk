#ifndef _FTK_DISTRIBUTED_UNION_FIND_H
#define _FTK_DISTRIBUTED_UNION_FIND_H

#include <vector>
#include <map>
#include <set>
#include <iostream>

// Union-find algorithm with distributed-memory parallelism
  // All elements need to be added by invoking the add function at the initialization stage. 

// Reference 
  // Paper: "Evaluation of connected-component labeling algorithms for distributed-memory systems"

// Add the sparse representation by using Hash map/table 

namespace ftk {
template <class IdType=std::string>
struct distributed_union_find
{
  distributed_union_find() {
    
  }

  // Initialization
  // Add and initialize elements
  void add(IdType i) {
    eles.insert(i); 
    id2parent[i] = i;
    sz[i] = 1;
  }

  // Operations

  void set_parent(IdType i, IdType par) {
    if(!has(i)) {
      throw "No such element. ";
    }

    id2parent[i] = par; 
  }


  // Queries

  bool has(IdType i) {
    return eles.find(i) != eles.end(); 
  }

  IdType parent(IdType i) {
    if(!has(i)) {
      throw "No such element. ";
    }

    return id2parent[i]; 
  }

  bool is_root(IdType i) {
    if(!has(i)) {
      throw "No such element. ";
    }

    return i == id2parent[i]; 
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
