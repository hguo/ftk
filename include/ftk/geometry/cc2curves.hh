#ifndef _FTK_CC2CURVE_HH
#define _FTK_CC2CURVE_HH

#include <ftk/algorithms/cca.hh>
#include <list>
#include <set>

namespace ftk {

template <typename NodeType>
std::vector<std::vector<NodeType>> connected_component_to_linear_components(
    const std::set<NodeType>& connected_component,
    const std::function<std::set<NodeType>(NodeType)>& neighbors)
{
  std::vector<std::vector<NodeType>> linear_components;
  if (connected_component.empty())
    return linear_components;

  std::set<NodeType> ordinary_nodes, special_nodes;
  for (const auto node : connected_component) {
    std::vector<NodeType> valid_neighbors;
    for (const auto neighbor : neighbors(node))
      if (neighbor != node && connected_component.find(neighbor) != connected_component.end())
        valid_neighbors.push_back(neighbor);

    if (valid_neighbors.size() > 2) {
      // std::cerr << "deg=" << valid_neighbors.size() << "," << node << std::endl;
      special_nodes.insert(node);
    }
    else ordinary_nodes.insert(node);
  }

  // connected components of ordinary nodes
  auto cc = extract_connected_components<NodeType, std::set<NodeType>>(
      [connected_component, neighbors, special_nodes](NodeType node) { // neighbor function
        auto my_neighbors = neighbors(node);
        for (auto it = my_neighbors.begin(); it != my_neighbors.end(); ) 
          if (connected_component.find(*it) == connected_component.end() 
              || special_nodes.find(*it) != special_nodes.end())
            it = my_neighbors.erase(it);
          else it ++;
        return my_neighbors;
      }, ordinary_nodes);

  // sort the linear graphs
  // fprintf(stderr, "#nodes=%lu, #linear graphs=%lu\n", connected_component.size(), cc.size());
  for (auto &c : cc) {
    std::list<NodeType> trace;
    std::set<NodeType> visited;

    auto seed = *c.begin();
    visited.insert(seed);
    trace.push_back(seed);
    c.erase(c.begin());

    // fprintf(stderr, "#node=%lu\n", c.size());
    // std::cerr << "seed: " << seed << std::endl;

    std::set<NodeType> seed_neighbors; 
    for (auto my_neighbor : neighbors(seed))
      if (ordinary_nodes.find(my_neighbor) != ordinary_nodes.end())
        seed_neighbors.insert(my_neighbor);

    if (seed_neighbors.size() == 0) break;
    for (int dir = 0; dir < 2; dir ++) {
      NodeType current = dir == 0 ? (*seed_neighbors.begin()) : (*seed_neighbors.rbegin());
      bool first_iteration = true;
      // fprintf(stderr, "dir=%d\n", dir);
      while (1) {
        if (first_iteration) first_iteration = false;
        else {
          if (dir == 0) trace.push_back(current);
          else (trace.push_front(current));
        }

        visited.insert(current);
        c.erase(current);

        bool found_next = false;
        for (auto my_neighbor : neighbors(current)) {
          if (my_neighbor != current && 
              c.find(my_neighbor) != c.end() && 
              // connected_component.find(my_neighbor) != connected_component.end() &&
              // special_nodes.find(my_neighbor) != special_nodes.end() && // TODO: handle special nodes
              visited.find(my_neighbor) == visited.end()) 
          {
            found_next = true;
            current = my_neighbor;
            break;
          }
        }
        if (!found_next) break;
      }
      if (seed_neighbors.size() == 1) break; // only one direction available
    }
    // fprintf(stderr, "trace.size=%lu\n", trace.size());

    linear_components.push_back(std::vector<NodeType>({trace.begin(), trace.end()}));
  }
  
  return linear_components;
}

template <typename NodeType>
bool is_loop(const std::vector<NodeType>& linear_graph, std::function<std::set<NodeType>(NodeType)> neighbors)
{
  if (linear_graph.size() == 0) return false;
  else if (linear_graph.size() == 1) return true;
  else {
    const auto front_neighbors = neighbors(linear_graph.front());
    return front_neighbors.find(linear_graph.back()) != front_neighbors.end();
  }
}

}

#endif
