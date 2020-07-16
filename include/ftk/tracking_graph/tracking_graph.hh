#ifndef _FTK_TRACKING_GRAPH_HH
#define _FTK_TRACKING_GRAPH_HH

#include <map>
#include <set>
#include <unordered_set>
#include <queue>
#include <mutex>
#include <functional>
#include <fstream>
#include <ftk/algorithms/cca.hh>
#include <ftk/tracking_graph/event.hh>

namespace ftk {

template <class TimeIndexType=size_t, class LabelIdType=size_t, class GlobalLabelIdType=size_t, class WeightType=int>
class tracking_graph {
public:
  tracking_graph();

  std::vector<TimeIndexType> get_timesteps() const;

  bool has_global_label(TimeIndexType t, LabelIdType l) const;
  GlobalLabelIdType get_global_label(TimeIndexType t, LabelIdType l) const;

  bool has_node(TimeIndexType t, LabelIdType l) const;
  bool has_edge(TimeIndexType t0, LabelIdType l0, TimeIndexType t1, LabelIdType l1) const;
  
  void add_node(TimeIndexType t, LabelIdType l);
  void add_edge(TimeIndexType t0, LabelIdType l0, TimeIndexType t1, LabelIdType l1);

  const std::map<TimeIndexType, std::vector<Event<TimeIndexType, LabelIdType> > > &get_events() const {return events;}

  void detect_events();
  void detect_events(TimeIndexType t0, TimeIndexType t1);

  void relabel();

  void generate_dot_file(const std::string& filename) const;

private:
  typedef std::pair<TimeIndexType, LabelIdType> Node;
  
  std::map<TimeIndexType, std::set<Node> > nodes;
  std::map<Node, std::set<Node> > left_links, right_links;

  std::map<Node, GlobalLabelIdType> nodeToGlobalLabelMap;
  std::map<GlobalLabelIdType, std::set<Node> > globalLabelToNodeMap;

  std::map<TimeIndexType, std::vector<Event<TimeIndexType, LabelIdType> > > events;

private: // counters for global label
  std::function<void()> resetGlobalLabelCallback;
  std::function<GlobalLabelIdType()> newGlobalLabelCallback;

  size_t defaultGlobalLabelCounter = 0;
  void defaultResetGlobalLabels() {defaultGlobalLabelCounter = 0;}
  GlobalLabelIdType defaultNewGlobalLabel() {return GlobalLabelIdType(++defaultGlobalLabelCounter);}

private:
  mutable std::mutex mutex;
};

////////////////////////////////////////////




////////////////////////////////////////////
  
template <class TimeIndexType, class LabelIdType, class GlobalLabelIdType, class WeightType>
tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType>::tracking_graph() : 
  resetGlobalLabelCallback(std::bind(&tracking_graph::defaultResetGlobalLabels, this)),
  newGlobalLabelCallback(std::bind(&tracking_graph::defaultNewGlobalLabel, this))
{}

template <class TimeIndexType, class LabelIdType, class GlobalLabelIdType, class WeightType>
bool tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType>::has_global_label(TimeIndexType t, LabelIdType l) const 
{
  return nodeToGlobalLabelMap.find(std::make_pair(t, l)) != nodeToGlobalLabelMap.end();
}

template <class TimeIndexType, class LabelIdType, class GlobalLabelIdType, class WeightType>
GlobalLabelIdType tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType>::get_global_label(TimeIndexType t, LabelIdType l) const 
{
  auto it = nodeToGlobalLabelMap.find(std::make_pair(t, l));
  if (it != nodeToGlobalLabelMap.end()) 
    return it->second;
  else
    return GlobalLabelIdType(-1); // TODO
}

template <class TimeIndexType, class LabelIdType, class GlobalLabelIdType, class WeightType>
bool tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType>::has_node(TimeIndexType t, LabelIdType l) const 
{
  std::unique_lock<std::mutex> lock(mutex);
  return nodes.find(std::make_pair(t, l)) != nodes.end();
}
  
template <class TimeIndexType, class LabelIdType, class GlobalLabelIdType, class WeightType>
bool tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType>::has_edge(TimeIndexType t0, LabelIdType l0, TimeIndexType t1, LabelIdType l1) const 
{
  std::unique_lock<std::mutex> lock(mutex);
  
  auto it = right_links.find(std::make_pair(t0, l0)); 
  if (it == right_links.end()) 
    return false;
  else if (it->find(std::make_pair(t1, l1)) == it->end())
    return false;
  else return true;
}
  
template <class TimeIndexType, class LabelIdType, class GlobalLabelIdType, class WeightType>
void tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType>::add_node(TimeIndexType t, LabelIdType l) 
{
  std::unique_lock<std::mutex> lock(mutex);
  nodes[t].insert(std::make_pair(t, l));
}
  
template <class TimeIndexType, class LabelIdType, class GlobalLabelIdType, class WeightType>
void tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType>::add_edge(TimeIndexType t0, LabelIdType l0, TimeIndexType t1, LabelIdType l1) 
{
  std::unique_lock<std::mutex> lock(mutex);

  auto n0 = std::make_pair(t0, l0); 
  auto n1 = std::make_pair(t1, l1);

  nodes[t0].insert(n0);
  nodes[t1].insert(n1); 

  left_links[n1].insert(n0);
  right_links[n0].insert(n1);
}
  
template <class TimeIndexType, class LabelIdType, class GlobalLabelIdType, class WeightType>
void tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType>::relabel()
{
  std::set<Node> allNodes;
  for (const auto &kv : nodes) {
    allNodes.insert(kv.second.begin(), kv.second.end());
  }

  auto components = extract_connected_components<Node, std::set<Node> >(
      [this](Node n) {
        std::set<Node> neighbors;
        if (left_links[n].size() == 1) // nodes need to be simply connected
          neighbors.insert(left_links[n].begin(), left_links[n].end());
        if (right_links[n].size() == 1)
          neighbors.insert(right_links[n].begin(), right_links[n].end());
        return neighbors;
      }, allNodes);

  // fprintf(stderr, "#components=%zu\n", components.size());

  nodeToGlobalLabelMap.clear();
  globalLabelToNodeMap.clear();

  resetGlobalLabelCallback();
  for (auto component : components) {
    GlobalLabelIdType globalLabel = newGlobalLabelCallback();
    globalLabelToNodeMap[globalLabel] = component;
    for (auto node : component) {
      nodeToGlobalLabelMap[node] = globalLabel;
    }
  }
}

template <class TimeIndexType, class LabelIdType, class GlobalLabelIdType, class WeightType>
std::vector<TimeIndexType> tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType>::get_timesteps() const
{
  std::vector<TimeIndexType> timesteps;
  for (auto kv : nodes)
    timesteps.push_back(kv.first);

  return timesteps;
}

template <class TimeIndexType, class LabelIdType, class GlobalLabelIdType, class WeightType>
void tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType>::detect_events(TimeIndexType t0, TimeIndexType t1)
{
  std::set<Node> intervalNodes;
  for (auto n : nodes[t0]) intervalNodes.insert(n);
  for (auto n : nodes[t1]) intervalNodes.insert(n);

  auto intervalNeighbors = [this, t0](Node n) {
    if (n.first == t0) return right_links[n];
    else return left_links[n];
  };

  auto components = extract_connected_components<Node, std::set<Node> >(intervalNeighbors, intervalNodes);
  // fprintf(stderr, "====%zu\n", components.size());
  for (auto component : components) {
    if (component.size() != 2) {
      Event<TimeIndexType, LabelIdType> e;
      e.interval = std::make_pair(t0, t1);
      for (auto n : component) {
        if (n.first == t0) e.lhs.insert(n.second);
        else e.rhs.insert(n.second);
      }
      events[t0].push_back(e);
      
      // json j = e;
      // fprintf(stderr, "%s\n", j.dump().c_str());
    }
  }
}

template <class TimeIndexType, class LabelIdType, class GlobalLabelIdType, class WeightType>
void tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType>::detect_events()
{
  const auto timesteps = get_timesteps();
  if (timesteps.size() < 2) return; // no intervals available

  // iterator over intervals
  for (size_t i = 0; i < timesteps.size() - 1; i ++) {
    TimeIndexType t0 = timesteps[i], t1 = timesteps[i + 1];
    detect_events(t0, t1);
  }
}

template <class TimeIndexType, class LabelIdType, class GlobalLabelIdType, class WeightType>
void tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType>::generate_dot_file(const std::string& filename) const 
{
  using namespace std;
  ofstream ofs(filename.c_str());
  if (!ofs.is_open()) return;

  auto node2str = [this](Node n) {
    return "T" + std::to_string(n.first) + "L" + std::to_string(n.second) + "G" + std::to_string(get_global_label(n.first, n.second));
    // std::stringstream ss; 
    // ss << "T" << n.first << "L" << n.second << "G" << get_global_label(n.first, n.second);
    // return ss.str();
  };

  ofs << "digraph {" << endl;
  ofs << "ratio = compress;" << endl;
  ofs << "rankdir = LR;" << endl;
  ofs << "ranksep =\"1.0 equally\";" << endl;
  ofs << "node [shape=circle];" << endl;
  // ofs << "node [shape=point,width=0,height=0];" << endl;
 
  for (const auto &kv : nodes) {
    for (const auto &n : kv.second) {
      auto globalLabel = get_global_label(n.first, n.second);
      int c = globalLabel % 6;
      std::string color;
        
      if (c == 0) color = "blue";
      else if (c == 1) color = "green";
      else if (c == 2) color = "cyan";
      else if (c == 3) color = "red";
      else if (c == 4) color = "purple";
      else if (c == 5) color = "yellow";

      ofs << node2str(n)
          << " [style=filled, fillcolor=" << color << "];" << endl;
    }

    ofs << "{rank=same; ";
    for (auto n : kv.second) {
      ofs << node2str(n) << ",";
    }
    ofs.seekp(-1, std::ios_base::end);
    ofs << "}" << endl;
  }

  for (const auto &kv : right_links) {
    const auto lNode = kv.first;
    // const int weight = kv.second.size() == 1 ? TODO: colsum & rowsum
    for (const auto &rNode : kv.second) {
      ofs << node2str(lNode) << "->" 
          << node2str(rNode) << endl;
          // << " [weight = " << weight << "];" << endl;
    }
  }

#if 0
  // node colors
  for (const auto& kv : tr.labels()) {
  }
#endif 
  ofs << "}" << endl;
  ofs.close();
}

}
#endif
