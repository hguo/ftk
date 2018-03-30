#ifndef _FTK_TRACKING_GRAPH
#define _FTK_TRACKING_GRAPH

#include "ftk/base/object.h"

typedef std::pair<int, int> ftkNodeId; // <timestep, localId>
// typedef std::pair<ftkNodeId, ftkNodeId> ftkEdgeId;

struct ftkNodeInfo {
  int globalId;
  std::set<ftkNodeId> forwardNodeIds, backwardNodeIds;
};

class ftkTrackingGraph : public ftkObject {
public:
  // inline ftkNodeInfo& addNode(int t, int lid) {return graph[makeNode(t, lid)];};
  inline void addNode(int t, int lid) {
    std::unique_lock<std::mutex> mlock(mutex);
    graph[makeNodeId(t, lid)];
  };
    
  inline void addEdge(int t0, int lid0, int t1, int lid1) {
    if (t0 == t1) return;
    
    std::unique_lock<std::mutex> mlock(mutex);
    ftkNodeId n0 = makeNodeId(t0, lid0), n1 = makeNodeId(t1, lid1);
    ftkNodeInfo &nodeInfo0 = graph[n0], &nodeInfo1 = graph[n1];

    if (t0 < t1) {
      nodeInfo0.forwardNodeIds.insert(n1);
      nodeInfo1.backwardNodeIds.insert(n0);
    } else {
      nodeInfo0.backwardNodeIds.insert(n1);
      nodeInfo1.forwardNodeIds.insert(n0);
    }
  }

public: // TODO: json
  inline json toJson() const {
    json j;
    return j;
  }

  inline void fromJson(json j) {

  }

private:
  std::map<ftkNodeId, ftkNodeInfo> graph;
  std::mutex mutex;
 
private:
  inline static ftkNodeId makeNodeId(int t, int lid) {return std::make_pair(t, lid);} 
};

//// 

#endif
