#ifndef _FTK_CC2CURVE_HH
#define _FTK_CC2CURVE_HH

#include <ftk/algorithms/cca.hh>
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
      if (connected_component.find(neighbor) != connected_component.end())
        valid_neighbors.push_back(neighbor);

    if (valid_neighbors.size() > 2) special_nodes.insert(node);
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
  for (auto &c : cc) {
    std::list<NodeType> trace;
    std::set<NodeType> visited;

    auto seed = *c.begin();
    visited.insert(seed);
    trace.push_back(seed);

    std::cerr << "seed: " << seed << std::endl;

    const auto seed_neighbors = neighbors(seed);
    if (seed_neighbors.size() == 0) break;

    for (int dir = 0; dir < 2; dir ++) {
      NodeType current = dir == 0 ? (*seed_neighbors.begin()) : (*seed_neighbors.rbegin());
      fprintf(stderr, "dir=%d\n", dir);
      while (1) {
        std::cerr << "current: " << current << std::endl;
        
        if (dir == 0) trace.push_back(current);
        else (trace.push_front(current));
        visited.insert(current);
        c.erase(current);

        bool found_next = false;
        for (auto my_neighbor : neighbors(current)) {
          if (c.find(my_neighbor) != c.end() && 
              // special_nodes.find(my_neighbor) != special_nodes.end() && // TODO: handle special nodes
              visited.find(my_neighbor) == visited.end()) 
          {
            found_next = true;
            current = my_neighbor;
            break;
          }
          else continue;
        }
        if (!found_next) break;
      }
      if (seed_neighbors.size() == 1) break; // only one direction available
    }
    fprintf(stderr, "trace.size=%lu\n", trace.size());
  }
}

}

#endif
