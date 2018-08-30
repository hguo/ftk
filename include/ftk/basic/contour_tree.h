#ifndef _FTK_CONTOUR_TREE_H
#define _FTK_CONTOUR_TREE_H

#include <set>
#include <map>

namespace ftk {

template <class IdType>
struct contour_tree {
  void add_node(IdType i) { nodes.insert(i); }
  bool has_node(IdType i) const { return nodes.find(i) != nodes.end(); }

  void add_arc(IdType lo, IdType hi) {
    upper_links[lo].insert(hi);
    lower_links[hi].insert(lo);
  }
  
  bool reduce_node(IdType i) { // remove a node whose up-degree and down-degree are both less or equal to 1
    size_t ud = upper_degree(i), ld = lower_degree(i);
    if (ud > 1 || ld > 1) return false; // cannot do the node reduction
    else if (ud == 1 && ld == 0) {
      lower_links[upper_node(i)].erase(i);
    } else if (ud == 0 && ld == 1) {
      upper_links[lower_node(i)].erase(i);
    } else { // ud == 1 && ld == 1
      auto hi = upper_node(i), lo = lower_node(i);
      lower_links[hi].erase(i);
      upper_links[lo].erase(i);
      lower_links[hi].insert(lo);
      upper_links[lo].insert(hi);
    }
    nodes.erase(i);
    upper_links.erase(i);
    lower_links.erase(i);
    return true;
  }
 
  bool is_leaf(IdType i) const {
    return upper_degree(i) + lower_degree(i) == 1;
  }

  size_t upper_degree(IdType lo) const {
    auto it = upper_links.find(lo);
    if (it != upper_links.end()) return it->second.size();
    else return 0;
  }

  size_t lower_degree(IdType hi) const {
    auto it = lower_links.find(hi);
    if (it != lower_links.end()) return it->second.size();
    else return 0;
  }

  IdType upper_node(IdType lo) const { // return the first upper node
    auto it = upper_links.find(lo);
    if (it != upper_links.end()) return *it->second.begin();
    else return IdType(-1);
  }

  IdType lower_node(IdType hi) const { // return the first lower node
    auto it = lower_links.find(hi);
    if (it != lower_links.end()) return *it->second.begin();
    else return IdType(-1);
  }

  const std::set<IdType>& upper_nodes(IdType i) const {
    auto it = upper_links.find(i);
    if (it != upper_links.end()) return it.second;
    else return std::set<IdType>();
  }

  const std::set<IdType>& lower_nodes(IdType i) const {
    auto it = lower_links.find(i);
    if (it != lower_links.end()) return it.second;
    else return std::set<IdType>();
  }

  void print() const {
    fprintf(stdout, "graph {\n");
    for (const auto &kv : upper_links) {
      IdType i = kv.first;
      for (const auto j : kv.second) 
        fprintf(stdout, "%zu -- %zu\n", i, j);
    }
    fprintf(stdout, "}\n");
  }

  template <class ValueType>
  void print_with_values(const std::function<ValueType(IdType)>& value) const {
    for (const auto &kv : upper_links) {
      IdType i = kv.first;
      for (const auto j : kv.second) {
        fprintf(stderr, "%zu(%f) -- %zu(%f);\n", 
            i, value(i), j, value(j));
      }
    }
  }

  void reduce() { // remove all nodes whose up-degree and down-degree are both 1
    std::vector<IdType> to_remove;
    for (auto i : nodes) 
      if (upper_degree(i) == 1 && lower_degree(i) == 1) {
        to_remove.push_back(i);
      }

    for (auto i : to_remove) {
      reduce_node(i);
    }
  }

protected:
  std::set<IdType> nodes;
  std::map<IdType, std::set<IdType> > upper_links, lower_links;
};

} // namespace ftk

#endif
