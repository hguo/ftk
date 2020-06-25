#ifndef _FTK_CONNECT_COMPONENT_TRACKER
#define _FTK_CONNECT_COMPONENT_TRACKER

#include <ftk/filters/filter.hh>
#include <ftk/tracking_graph/tracking_graph.hh>

namespace ftk {

struct connected_component_tracker : public filter
{
  connected_component_tracker() {}
  virtual ~connected_component_tracker() {};

  virtual bool advance_timestep();
  virtual void update_timestep();
  void update() {};
  void finalize();

  void push_labeled_data_snapshot(const std::vector<int>& labels);
  
  bool pop_snapshot();

  const ftk::tracking_graph<>& get_tracking_graph() const {return tg;}

protected:
  ftk::tracking_graph<> tg;
  std::deque<std::vector<int>> labeled_data_snapshots;
  int current_timestep = 0;
};

/////

inline void connected_component_tracker::finalize()
{
  tg.relabel();
}

inline void connected_component_tracker::update_timestep()
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

inline bool connected_component_tracker::advance_timestep()
{
  // fprintf(stderr, "#labeled_data_snapshots=%zu\n", labeled_data_snapshots.size());
  update_timestep();

  current_timestep ++;
  
  if (labeled_data_snapshots.size() > 1)
    pop_snapshot();
  return labeled_data_snapshots.size() > 0;
}

inline void connected_component_tracker::push_labeled_data_snapshot(const std::vector<int>& labels)
{
  labeled_data_snapshots.push_back(labels);
}

inline bool connected_component_tracker::pop_snapshot()
{
  if (labeled_data_snapshots.size() > 0) {
    labeled_data_snapshots.pop_front();
    return true;
  } else return false;
}

}

#endif
