#ifndef _MESH_GRAPH_HH
#define _MESH_GRAPH_HH

#include <set>
#include <vector>

namespace ftk {

template <typename index_type=size_t, typename chirality_type=signed char>
struct mesh_graph {
  virtual index_type n(int d) const = 0; // number of d-elements
  virtual bool valid(int d, index_type i) = 0; // check if d-element i is valid

  // children of mesh elements
  virtual std::vector<std::pair<index_type, chirality_type> > links_edge_node(index_type i) = 0; // returns a list of nodes with chirality=1
  virtual std::vector<std::pair<index_type, chirality_type> > links_face_edge(index_type i) = 0; // returns a list of edges with chirality=1
  virtual std::vector<std::pair<index_type, chirality_type> > links_cell_face(index_type i) = 0; // returns a list of faces with chirality=1
  
  virtual std::vector<std::pair<index_type, chirality_type> > children(int d, index_type i) { // d-1
    if (d == 1) return links_edge_node(i); 
    else if (d == 2) return links_face_edge(i);
    else if (d == 3) return links_cell_face(i);
    else return std::vector<std::pair<index_type, chirality_type>>();
  }

  // parents of mesh elements
  virtual std::vector<std::pair<index_type, chirality_type> > links_node_edge(index_type i) = 0;
  virtual std::vector<std::pair<index_type, chirality_type> > links_edge_face(index_type i) = 0;
  virtual std::vector<std::pair<index_type, chirality_type> > links_face_cell(index_type i) = 0;
  
  virtual std::vector<std::pair<index_type, chirality_type> > parents(int d, index_type i) { // d+1
    if (d == 0) return links_node_edge(i);
    else if (d == 1) return links_edge_face(i);
    else if (d == 2) return links_face_cell(i);
    else return std::vector<std::pair<index_type, chirality_type>>();
  }


  // derived links that are optional
  virtual std::vector<std::pair<index_type, chirality_type> > links_cell_nodes(index_type i) {
    return std::vector<std::pair<index_type, chirality_type> >();
  }

  virtual std::vector<std::pair<index_type, chirality_type> > links_cell_edges(index_type i) { // usually not useful
    return std::vector<std::pair<index_type, chirality_type> >();
  }

  virtual std::vector<std::pair<index_type, chirality_type>> links_face_nodes(index_type i) {
    std::vector<std::pair<index_type, chirality_type>> results;
    const auto edges = links_face_edge(i);
    for (auto i = 0; i < edges.size(); i ++) {
      const auto edge = edges[i];
      const auto nodes = links_edge_node(edge.first);
      results.push_back(edge.second == 1 ? nodes[1] : nodes[0]);
    }
    return results;
  }
};

}

#endif
