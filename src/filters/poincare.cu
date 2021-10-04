#include "common.cuh"
#include "mx4.cuh"
#include <ftk/numeric/rk4.hh>

using namespace ftk;

template <typename I, typename F>
__device__
bool 
