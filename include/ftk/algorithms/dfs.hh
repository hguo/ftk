#ifndef _DFS_H
#define _DFS_H

#include <stack>
#include <set>
#include <functional>

namespace ftk {

template <class Node, class ContainerType>
void dfs(
    Node seed, 
    const std::function<ContainerType(Node)> &neighbors,
    const std::function<void(Node)> &operation,
    std::function<bool(Node)> criterion = [](Node){return true;}) 
{
  if (!criterion(seed)) return;

  std::set<Node> visited;
  std::stack<Node> S;
  S.push(seed);
  visited.insert(seed);

  while (!S.empty()) {
    Node current = Q.front(); // current = S.pop();
    Q.pop();

    operation(current);

    for (auto n : neighbors(current)) {
      if (visited.find(n) == visited.end()) // not visited
        if (criterion(n)) {
          visited.insert(n);
          Q.push(n); // S.push(n);
        }
    }
  }
}

}

#endif
