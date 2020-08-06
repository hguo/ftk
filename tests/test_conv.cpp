#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <vector>
#include <cmath>
#include <ftk/ndarray/conv.hh>
#include <ftk/ndarray.hh>
#include <ftk/numeric/rand.hh>

const int nruns = 100000;
const double epsilon = 1e-9;

TEST_CASE("2D_conv_test") {
  ftk::ndarray<double> k({2,2}),d({3,3}),r;
  std::vector<double> ans = {17,21,29,33,},
                      res;
  k.from_vector({-1,2,1,2,});
  d.from_vector({1,2,3,4,5,6,7,8,9,});
  r = ftk::conv2D(d,k);
  r.to_vector(res);
  for (int i=0;i<4;++i){
    ans[i]/=4;
  }
  REQUIRE(res == ans);
}

TEST_CASE("2D_conv_padding_test") {
  ftk::ndarray<double> k({3,3}),d({3,3}),r;
  std::vector<double> ans = {-13,-20,-17,-18,-24,-18,13,20,17,},
                      res;
  k.from_vector({1,2,1,0,0,0,-1,-2,-1,});
  d.from_vector({1,2,3,4,5,6,7,8,9,});
  r = ftk::conv2D(d,k,1);
  r.to_vector(res);
  for (int i=0;i<9;++i){
    ans[i]/=9;
  }
  REQUIRE(res == ans);
}

TEST_CASE("2D_conv_gaussian_test") {
  ftk::ndarray<double> d({3,3}),r;
  std::vector<double> ans={3/4.,1.,6/4.,7/4.,},
                      res;
  d.from_vector({1,2,3,4,5,6,7,8,9,});
  r = ftk::conv2D_gaussian(d,1.,2,2);
  r.to_vector(res);
  REQUIRE(res == ans);
}

TEST_CASE("3D_conv_test") {
  ftk::ndarray<double> d({3,3,3}),k({2,2,2}),r;
  std::vector<double> ans = {-4, -3.875, -3.625, -3.5, -2.875, -2.75, -2.5, -2.375,},
                      data,res;
  for (int i=0;i<27;i++)
    data.push_back(i+1);
  k.from_vector({1,2,1,0,0,0,-1,-2,});
  d.from_vector(data);
  r = ftk::conv3D(d,k,0);
  r.to_vector(res);
  REQUIRE(res == ans);
}

TEST_CASE("3D_conv_padding_test") {
  ftk::ndarray<double> d({2,2,2}),k({2,2,2}),r;
  std::vector<double> ans = { -0.25,-0.625,-0.25,-0.75,-1.375,-0.5,0,0,0,-1.25,-2,-0.5,-1.5,-1.875,-0.25,0.75,1.375,0.5,0,0.625,0.75,1.25,3,1.75,1.75,2.875,1,},
                      data,res;
  for (int i=0;i<8;i++)
    data.push_back(i+1);
  k.from_vector({1,2,1,0,0,0,-1,-2,});
  d.from_vector(data);
  r = ftk::conv3D(d,k,1);
  r.to_vector(res);
  REQUIRE(res == ans);
}

TEST_CASE("3D_conv_gaussian_test") {
  ftk::ndarray<double> d({3,3,3}),r({2,2,2});
  std::vector<double> ans={-0.6875, -0.8125, -1.0625, -1.1875, -1.8125, -1.9375, -2.1875, -2.3125},
                      data,
                      res;
  for (int i=0;i<27;i++)
    data.push_back(-i+1);
  d.from_vector(data);
  r = ftk::conv3D_gaussian(d,1.,2,2,2);
  r.to_vector(res);
  REQUIRE(res == ans);
}

int main(int argc, char **argv)
{
  Catch::Session session;
  return session.run(argc, argv);
}
