#ifndef _FTK_DIJKSTRA_H
#define _FTK_DIJKSTRA_H

#include <vector>
#include <set>
#include <map>
#include <limits>
#include <algorithm>
#include <functional>

namespace ftk {

template <typename Node, typename T>
std::vector<Node> 
dijkstra(
    const Node src, const Node dst,
    const std::function<std::set<Node>(Node)> neighbor_func,
    const std::function<T(Node, Node)> dist_func)
{
  std::set<Node> visited;
  std::map<Node, T> dist;
  std::map<Node, Node> prev;

  dist[src] = 0;

  while (1) {
    Node u; // vertex in Q w/ min dist[u]
    T mindist = std::numeric_limits<T>::max();
    for (const auto &kv : dist)
      if (visited.find(kv.first) == visited.end() && kv.second < mindist) {
        mindist = kv.second;
        u = kv.first;
      }

    visited.insert(u);
    // fprintf(stderr, "erasing %d, mindist=%f, #Q=%zu\n", u, mindist, Q.size());
    if (u == dst) break;

    T alt;
    for (const auto v : neighbor_func(u)) {
      alt = dist[u] + dist_func(u, v);
      bool change = false;
      if (dist.find(v) == dist.end()) // dist[v]=inf
        change = true;
      else if (alt < dist[v])
        change = true;

      if (change) {
        dist[v] = alt;
        prev[v] = u;
      }
    }
  }

  // building sequence
  std::vector<Node> result;
  Node u = dst;
  while (1) {
    result.push_back(u);
    
    if (u == src || prev.find(u) == prev.end()) break;
    else u = prev[u];
  }
  std::reverse(result.begin(), result.end());

  return result;
}

}

#endif
