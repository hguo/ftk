#ifndef _FTK_CONNECT_COMPONENT_TRACKER
#define _FTK_CONNECT_COMPONENT_TRACKER

#include <ftk/filters/filter.hh>
#include <ftk/tracking_graph/tracking_graph.hh>

namespace ftk {

template <typename TimeIndexType=size_t, typename LabelIdType=size_t>
struct connected_component_tracker : public filter
{
  connected_component_tracker() {}
  virtual ~connected_component_tracker() {};

  virtual bool advance_timestep();
  virtual void update_timestep();
  void update() {};
  void finalize();

  virtual void push_labeled_data_snapshot(const std::vector<LabelIdType>& labels);
  const std::vector<LabelIdType>& get_last_labeled_data_snapshot() const {return labeled_data_snapshots.back();}
  
  template <typename ContainerType>
  void push_unlabeled_data_snapshot(const std::vector<LabelIdType>& labels, std::function<ContainerType(LabelIdType)> neighbors);

  bool pop_snapshot();

  const ftk::tracking_graph<>& get_tracking_graph() const {return tg;}

protected:
  ftk::tracking_graph<TimeIndexType, LabelIdType> tg;
  std::deque<std::vector<LabelIdType>> labeled_data_snapshots;
  TimeIndexType current_timestep = 0;
};

/////

template <typename TimeIndexType, typename LabelIdType>
void connected_component_tracker<TimeIndexType, LabelIdType>::finalize()
{
  tg.relabel();
}

template <typename TimeIndexType, typename LabelIdType>
void connected_component_tracker<TimeIndexType, LabelIdType>::update_timestep()
{
  if (labeled_data_snapshots.size() < 2) return;

  const auto &labels0 = labeled_data_snapshots[0],
             &labels1 = labeled_data_snapshots[1];

  for (size_t i = 0; i < labels0.size(); i ++)
    if (labels0[i] != 0 && labels1[i] != 0) {
      // fprintf(stderr, "edge: %d --> %d\n", labels0[i], labels1[i]);
      tg.add_edge(current_timestep-1, labels0[i], current_timestep, labels1[i]);
    }
}

template <typename TimeIndexType, typename LabelIdType>
bool connected_component_tracker<TimeIndexType, LabelIdType>::advance_timestep()
{
  // fprintf(stderr, "#labeled_data_snapshots=%zu\n", labeled_data_snapshots.size());
  update_timestep();

  current_timestep ++;
  
  if (labeled_data_snapshots.size() > 1)
    pop_snapshot();
  return labeled_data_snapshots.size() > 0;
}

template <typename TimeIndexType, typename LabelIdType>
void connected_component_tracker<TimeIndexType, LabelIdType>::push_labeled_data_snapshot(const std::vector<LabelIdType>& labels)
{
  labeled_data_snapshots.push_back(labels);
}

template <typename TimeIndexType, typename LabelIdType>
bool connected_component_tracker<TimeIndexType, LabelIdType>::pop_snapshot()
{
  if (labeled_data_snapshots.size() > 0) {
    labeled_data_snapshots.pop_front();
    return true;
  } else return false;
}

}

#endif
