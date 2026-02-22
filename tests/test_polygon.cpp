#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <ftk/numeric/sign_det.hh>

const int n1 = 7;
const int indices1[n1] = {1, 2, 3, 4, 5, 6, 7};
const int verts1[n1][2] = {
  {10000, 10000},
  {0, 10000}, 
  {-10000, 0},
  {-10000, -10000},
  {0, -10000},
  {10000, -10000},
  {10000, 0}
};

TEST_CASE("point_in_polygon2") {
  const int x1[2] = {1, 1};
  REQUIRE(ftk::robust_point_in_polygon2<int>(
        x1, n1, indices1, verts1) == true);
  
  const int x2[2] = {20000, 0};
  REQUIRE(ftk::robust_point_in_polygon2<int>(
        x2, n1, indices1, verts1) == false);
}

#include "main.hh"
