#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include <ftk/algorithms/hoshen_kopelman.hh>
#include <random>

TEST_CASE("hoshen_kopelman_2d", "hoshen_kopelman")
{
  ftk::ndarray<int> input({8, 8}), output({8, 8});
  input.fill({
    1,1,1,1,1,1,1,1,
    0,0,0,0,0,0,0,1,
    1,0,0,0,0,1,0,1,
    1,0,0,1,0,1,0,1,
    1,0,0,1,0,1,0,1,
    1,0,0,1,1,1,0,1,
    1,1,1,1,0,0,0,1,
    0,0,0,1,1,1,0,1});
  output.fill({
    1,1,1,1,1,1,1,1, 
    0,0,0,0,0,0,0,1, 
    2,0,0,0,0,2,0,1, 
    2,0,0,2,0,2,0,1, 
    2,0,0,2,0,2,0,1, 
    2,0,0,2,2,2,0,1, 
    2,2,2,2,0,0,0,1, 
    0,0,0,2,2,2,0,1});

  SECTION("test hoshen kopelman 2d") {
    int nc = hoshen_kopelman_2d(input);
    CHECK(nc == 2);
    CHECK(input == output);
  }
}

#include "main.hh"
