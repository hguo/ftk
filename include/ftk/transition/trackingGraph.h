#ifndef _FTK_TRACKING_GRAPH_H
#define _FTK_TRACKING_GRAPH_H

#include "ftk/base/object.h"
#include <utility>
#include <map>
#include <set>
#include <queue>
#include <mutex>

struct ftkNode {
  int gid = -1; // global id.  gid is -1 if not labeled
  std::pair<int, int> id; // <time, lid>
  std::map<int, std::set<int> > neighbors; // <time, lids>
};

struct ftkTrackingGraph { 
  inline void addNode(int t, int lid);
  inline void addEdge(int t0, int lid0, int t1, int lid1);

  inline void relabel();
  inline void detectEvents();

  std::map<int, std::map<int, ftkNode> > nodes; // <time, <lid, node> >
  std::set<std::pair<int, int> > unlabeledNodeIds; // <time, lid>
  std::map<int, std::map<int, int> > labeledComponents; // <gid, <time, lid> >
  std::mutex mutex;
};


///////////////////

void ftkTrackingGraph::detectEvents()
{
  if (!unlabeledNodeIds.empty()) relabel(); 

  // extract endpoints of components
  std::set<std::pair<int, int> > endpoints; // <time, lid>
  for (const auto &kv : labeledComponents) {
    endpoints.insert(*kv.second.begin());
    endpoints.insert(*kv.second.rbegin());
  }

  // detect connected comonents in endpoints
  std::set<std::set<std::pair<int, int> > > connectedEndpoints;
  std::set<std::pair<int, int> > visitedEndpoints;
  std::queue<std::pair<int, int> > Q;

  Q.push(*endpoints.begin());
  while (!Q.empty()) {
    std::pair<int, int> currentNodeId = Q.front(); 
    const ftkNode& currentNode = nodes[currentNodeId.first][currentNodeId.second];
    visitedEndpoints.insert(currentNodeId);
    Q.pop();

    // find endpoint neighbors


  }

  for (const auto &nodeId : endpoints) {
    const ftkNode& node = nodes[nodeId.first][nodeId.second];
  }
}

void ftkTrackingGraph::relabel()
{
  int nComponents = 0;
  
  while (!unlabeledNodeIds.empty()) {
    // fprintf(stderr, "new component!\n");
    std::map<int, int> component;

    std::queue<std::pair<int, int> > Q;
    Q.push(*unlabeledNodeIds.begin());
    // Q.push(std::make_pair(nodes.begin()->first, nodes.begin()->second.begin()->first));

    while (!Q.empty()) {
      std::pair<int, int> currentNodeId = Q.front(); // time, lid
      const ftkNode& currentNode = nodes[currentNodeId.first][currentNodeId.second];
      Q.pop();

      // process current node
      fprintf(stderr, "%d, %d\n", currentNodeId.first, currentNodeId.second);
      component.insert(currentNodeId);
      unlabeledNodeIds.erase(currentNodeId);

      for (const auto &kv : currentNode.neighbors) {
        if (kv.second.size() == 1) { // simple connected component
          for (const auto &neighbor : kv.second) {
            std::pair<int, int> neighborId(kv.first, neighbor);
            if (unlabeledNodeIds.find(neighborId) != unlabeledNodeIds.end()) // neighbor is unlabeled
              Q.push(neighborId);
          }
        }
      }
    }

    labeledComponents[nComponents++] = component;
  }
}

void ftkTrackingGraph::addNode(int t, int lid) {
  std::unique_lock<std::mutex> mlock(mutex);
  ftkNode &node = nodes[t][lid];
  node.id = std::make_pair(t, lid);
  unlabeledNodeIds.insert(std::make_pair(t, lid));
}

void ftkTrackingGraph::addEdge(int t0, int lid0, int t1, int lid1)
{
  std::unique_lock<std::mutex> mlock(mutex);
  ftkNode &node0 = nodes[t0][lid0];
  node0.neighbors[t1].insert(lid1);

  ftkNode &node1 = nodes[t1][lid1];
  node1.neighbors[t0].insert(lid0);
}

#endif
