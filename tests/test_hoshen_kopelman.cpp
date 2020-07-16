#include <gtest/gtest.h>
#include <ftk/algorithms/hoshen_kopelman.hh>
#include <random>

class hoshen_kopelman_test : public testing::Test {
public:
};

TEST_F(hoshen_kopelman_test, hoshen_kopelman_2d)
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

  int nc = hoshen_kopelman_2d(input);
  EXPECT_EQ(nc, 2);
  EXPECT_EQ(input, output);
}
