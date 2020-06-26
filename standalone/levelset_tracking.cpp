#include <fstream>
#include <mutex>
#include <set>
#include <cassert>
#include <ftk/external/cxxopts.hpp>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/filters/levelset_tracker.hh>
#include <ftk/filters/connected_component_tracker.hh>
#include <ftk/ndarray.hh>
#include <ftk/algorithms/hoshen_kopelman.hh>
#include <ftk/tracking_graph/tracking_graph.hh>
#include "cli_constants.hh"

const size_t DW = 64, DH = 32, DT = 128;
const double threshold = 0.6;

ftk::ndarray<double> request_timestep(int k) // requesting k-th timestep
{
  // const double t = DT == 1 ? 0.0 : double(k)/(DT-1);
  // return ftk::synthetic_woven_2D<double>(DW, DH, t);
  const double t = DT == 1 ? 0.0 : double(k)/(DT-1) * 10;
  return ftk::synthetic_merger_2D<double>(DW, DH, t);
}

void track_levelset()
{
  auto *tracker = new ftk::levelset_tracker<>; // ftk::connected_component_tracker<>;
  tracker->set_threshold( threshold );

  for (int current_timestep = 0; current_timestep < DT; current_timestep ++) {
    fprintf(stderr, "current_timestep=%d\n", current_timestep);
    ftk::ndarray<double> field_data = request_timestep(current_timestep);
    // ftk::ndarray<size_t> label_data = threshold_filter<size_t>(field_data);
    // size_t nc = ftk::hoshen_kopelman_2d(label_data);
    // tracker->push_labeled_data_snapshot(label_data.std_vector());
    
    tracker->push_scalar_field_data_snapshot(field_data);
    tracker->advance_timestep();
  }

  tracker->finalize();

  const auto &tg = tracker->get_tracking_graph();
  tg.generate_dot_file("dot");

  delete tracker;
}

int main(int argc, char **argv)
{
  track_levelset();
  return 0;
}
