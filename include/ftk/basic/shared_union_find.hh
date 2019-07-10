#ifndef _FTK_SHARED_UNION_FIND_H
#define _FTK_SHARED_UNION_FIND_H

#include <vector>
#include <map>
#include <set>
#include <iostream>


// TBD
  // Complete thread-safe add function
  // Add unit test. 

// Union-find algorithm with shared-memory parallelism
  // All elements need to be added by invoking the add function at the initialization stage. 
  // This prograpm cannot be used for distributed-memory parallelism. 
  // This program assumes all elements are stored on the same memory. 

// Reference 
  // Paper: "A Randomized Concurrent Algorithm for Disjoint Set Union"

namespace ftk {
template <class IdType=std::string>
struct shared_union_find
{
  shared_union_find() : eles(), id2parent(), sz() {
    
  }

  // // Initialization
  // // Add and initialize elements
  // void add(IdType i) {
  //   eles.insert(i); 
  //   // id2parent[i] = i; 
  //   id2parent.insert(std::make_pair (i, i)); 
  //   sz[i] = 1;
  // }

  // Operations
  
  // Compare-and-swap
  	// Reference: https://en.wikipedia.org/wiki/Compare-and-swap
  bool CAS(std::pair<const IdType, IdType> ele_parent, IdType old_parent, IdType new_parent) {
    if(ele_parent->second != old_parent) {
      return false; 
    }

    ele_parent->second = new_parent; 

    return true; 
  }


  // Union by size
    // Return wehther i and j are in the same set. 
  bool unite(IdType i, IdType j) {
    // if(!has(i) || !has(j)) {
    //   throw "No such element. ";

    //   return false ;
    // }

    IdType u = i; 
    IdType v = j;

    while(true) {
      u = find(u); 
      v = find(v); 
      if(u < v) {
        if(CAS(id2parent.find(u), u, v)) {
          return false; 
        }
      } else if (u == v) {
        return true; 
      } else {
        if(CAS(id2parent.find(v), v, u)) {
          return false; 
        }
      }
    }

  }

  // // Queries

  bool has(IdType i) {
    return eles.find(i) != eles.end(); 
  }

  // IdType parent(IdType i) {
  //   if(!has(i)) {
  //     throw "No such element. ";
  //   }

  //   // return id2parent[i]; 
  //   return id2parent.find(i)->second; 
  // }

  // Find the root of an element. 
    // If the element is not in the data structure, return the element itself. 
    // Path compression by path halving method
  IdType find(IdType i) {
    // if(!has(i)) {
    //   throw "No such element. ";
    // }

    IdType u = i; 
    while(true) {
      IdType v = id2parent.find(u)->second; 
      IdType w = id2parent.find(v)->second; 
      if(v == w) {
      	return v; 
      } else {
      	CAS(id2parent.find(u), v, w); 
      	u = id2parent.find(u)->second;
      }
    }
  }

  // bool is_root(IdType i) {
  //   if(!has(i)) {
  //     throw "No such element. ";
  //   }

  //   return i == id2parent[i]; 
  // }

  bool same_set(IdType i, IdType j) {
  	if(!has(i) || !has(j)) {
      throw "No such element. ";
    }

  	IdType u = i; 
  	IdType v = j; 
  	while(true) {
  		u = find(u);
  		v = find(v); 

  		if(u == v) {
  			return true;
  		}

  		if(u == id2parent.find(u)->second) {
  			return false; 
  		}
  	}
  }

  // std::vector<std::set<IdType>> get_sets() {
  //   std::map<IdType, std::set<IdType>> root2set; 
  //   for(auto ite = eles.begin(); ite != eles.end(); ++ite) {
  //     IdType root = find(*ite); 

  //     root2set[root].insert(*ite);
  //   }

  //   std::vector<std::set<IdType>> results; 
  //   for(auto ite = root2set.begin(); ite != root2set.end(); ++ite) {
  //     // if(is_root(*ite))
  //     results.push_back(ite->second); 
  //   }

  //   return results; 
  // }

public:
  std::set<IdType> eles; 

private:
  // Use HashMap to support sparse union-find
  std::map<IdType, IdType> id2parent;
  std::map<IdType, size_t> sz;

};

}

#endif
