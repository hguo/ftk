#ifndef _FTK_GRAPH_H
#define _FTK_GRAPH_H

#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <ftk/basic/union_find.hh>

// Graph data structure based on union-find
  // Initialization
    // All vertices need to be added by using add_vertex
    // The vertices of an edge need to be added before adding the edge

namespace ftk {
template <class IdType=std::string>
struct graph
{
  graph() {
    
  }

  // Initialization
  // Add and initialize a vertex
  void add_vertex(IdType i) {
    uf.add(i); 
  }

  // Add and initialize an undirected edge
    // Assume each edge pair (i, j) is added once, i.e., we just need to add one of (i, j) and (j, i). 
      // Currently, do not remove duplicated vertices
  void add_edge(IdType i, IdType j) {
    uf.unite(i, j); 

    linked_verts_map[i].push_back(j); 
    linked_verts_map[j].push_back(i); 
  }

  // Queries

  bool has_vertex(IdType i) {
    return uf.has(i);
  }

  const std::vector<IdType> get_linked_vertices(IdType i) {
    return linked_verts_map[i]; 
  }

  std::vector<std::set<IdType>> get_connected_components() {
    return uf.get_sets(); 
  }

private:
  union_find<IdType> uf; 
  std::map<IdType, std::vector<IdType>> linked_verts_map;
};

}

#endif
