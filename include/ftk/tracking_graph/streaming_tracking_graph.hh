#ifndef _FTK_STREAMING_TRACKING_GRAPH_HH
#define _FTK_STREAMING_TRACKING_GRAPH_HH

#include <ftk/tracking_graph/tracking_graph.hh>

namespace ftk {

template <class TimeIndexType=size_t, class LabelIdType=size_t, class GlobalLabelIdType=size_t, class WeightType=int>
class streaming_tracking_graph : public tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType> {
public:
  streaming_tracking_graph();

  void add_node(LabelIdType); // add a node for the current timestep
  void add_edge(LabelIdType prev, LabelIdType curr); // add an edge between previous timestep and the current timestep
  void advance_timestep(); // advance timestep, update global labels, and detect events
};

////////////////////////////////////////////




////////////////////////////////////////////
template <class TimeIndexType, class LabelIdType, class GlobalLabelIdType, class WeightType>
void tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType>::add_node(LabelIdType l) const 
{
}

template <class TimeIndexType, class LabelIdType, class GlobalLabelIdType, class WeightType>
void tracking_graph<TimeIndexType, LabelIdType, GlobalLabelIdType, WeightType>::add_edge(LabelIdType l) const 
{
}

}
