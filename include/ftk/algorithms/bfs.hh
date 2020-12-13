#ifndef _BFS_H
#define _BFS_H

#include <queue>
#include <set>
#include <functional>

namespace ftk {

template <class Node, class ContainerType>
void bfs(
    Node seed, 
    const std::function<ContainerType(Node)> &neighbors,
    const std::function<void(Node)> &operation,
    std::function<bool(Node)> criterion = [](Node){return true;}) 
{
  if (!criterion(seed)) return;

  std::set<Node> visited;
  std::queue<Node> Q;
  Q.push(seed);
  visited.insert(seed);

  while (!Q.empty()) {
    Node current = Q.front();
    Q.pop();

    operation(current);

    for (auto n : neighbors(current)) {
      if (visited.find(n) == visited.end())
        if (criterion(n)) {
          visited.insert(n);
          Q.push(n);
        }
    }
  }
}

}

#endif
