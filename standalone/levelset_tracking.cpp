#include <fstream>
#include <mutex>
#include <set>
#include <cassert>
#include <ftk/external/cxxopts.hpp>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/filters/levelset_tracker_2d_regular.hh>
#include <ftk/ndarray.hh>
#include <ftk/algorithms/hoshen_kopelman.hh>
#include <ftk/tracking_graph/tracking_graph.hh>
#include "cli_constants.hh"

const size_t DW = 32, DH = 32, DT = 10;
const double threshold = 0.6;

ftk::ndarray<double> request_timestep(int k) // requesting k-th timestep
{
  const double t = DT == 1 ? 0.0 : double(k)/(DT-1);
  return ftk::synthetic_woven_2D<double>(DW, DH, t);
}

ftk::ndarray<int> threshold_filter(const ftk::ndarray<double>& array)
{
  ftk::ndarray<int> rtn;
  rtn.reshape(array);
  for (size_t i = 0; i < array.nelem(); i ++)
    rtn[i] = array[i] >= threshold ? 1 : 0;
  return rtn;
}

void track_levelset()
{
  ftk::levelset_tracker_2d_regular *tracker = new ftk::levelset_tracker_2d_regular;
  ftk::tracking_graph<> tg;

  int current_timestep = 0;
  ftk::ndarray<int> last_label_data;
  while (1) {
    ftk::ndarray<double> field_data = request_timestep(current_timestep);
    ftk::ndarray<int> label_data = threshold_filter(field_data);
    int nc = ftk::hoshen_kopelman_2d(label_data);
    for (int i = 0; i < nc; i ++)
      tg.add_node(current_timestep, i);

    if (current_timestep == DT - 1) { // last time step
      // tracker->update_timestep();
      break;
    } else if (current_timestep != 0) {
      for (size_t i = 0; i < label_data.nelem(); i ++) {
        if (last_label_data[i] && label_data[i])
          tg.add_edge(current_timestep-1, last_label_data[i], 
              current_timestep, label_data[i]);
      }
    }
    last_label_data = label_data;
    current_timestep ++;
  }

  tg.relabel();
  tg.generate_dot_file("dot");

  delete tracker;
}

int main(int argc, char **argv)
{
  track_levelset();
  return 0;
}
