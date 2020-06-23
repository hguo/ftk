#include <fstream>
#include <mutex>
#include <set>
#include <cassert>
#include <ftk/external/cxxopts.hpp>
#include <ftk/ndarray/synthetic.hh>
#include <ftk/filters/levelset_tracker_2d_regular.hh>
#include <ftk/ndarray.hh>
#include <ftk/algorithms/hoshen_kopelman.hh>
#include "cli_constants.hh"

const size_t DW = 32, DH = 32, DT = 10;

ftk::ndarray<double> request_timestep(int k) // requesting k-th timestep
{
  const double t = DT == 1 ? 0.0 : double(k)/(DT-1);
  return ftk::synthetic_woven_2D<double>(DW, DH, t);
}

void track_levelset()
{
  ftk::levelset_tracker_2d_regular *tracker = new ftk::levelset_tracker_2d_regular;
  delete tracker;
}

int main(int argc, char **argv)
{
  ftk::ndarray<int> input({8, 8});
  input.fill({
    1,1,1,1,1,1,1,1,
    0,0,0,0,0,0,0,1,
    1,0,0,0,0,1,0,1,
    1,0,0,1,0,1,0,1,
    1,0,0,1,0,1,0,1,
    1,0,0,1,1,1,0,1,
    1,1,1,1,0,0,0,1,
    0,0,0,1,1,1,0,1});

  std::cerr << input << std::endl;

  int nc = hoshen_kopelman(input);
  fprintf(stderr, "nc=%d\n", nc);

  input.print(std::cerr) << std::endl;

  return 0;
}
