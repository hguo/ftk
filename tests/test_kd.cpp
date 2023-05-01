#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <ftk/basic/kd.hh>

const int npts = 10000;
const int ntargets = 1000;

TEST_CASE("kd_nearest_neighbor2") {
  ftk::kd<double, 2> mykd;

  std::vector<double> pts(npts * 2);
  for (size_t i = 0; i < npts * 2; i ++)
    pts[i] = (double)rand() / RAND_MAX;

  mykd.set_inputs(npts, &pts[0]);
  mykd.build();

  for (auto batch = 0; batch < npts; batch ++) {
    std::array<double, 2> x = {(double)rand() / RAND_MAX, (double)rand() / RAND_MAX};

    size_t i = mykd.find_nearest(x);
    fprintf(stderr, "x=%f, %f\n", x[0], x[1]);
    fprintf(stderr, "%zu: %f, %f\n", i, mykd.pts[i][0], mykd.pts[i][1]);
    
    size_t j = mykd.find_nearest_naive(x);
    fprintf(stderr, "%zu: %f, %f\n", j, mykd.pts[j][0], mykd.pts[j][1]);

    REQUIRE(i == j);
  }
}

TEST_CASE("kd_nearest_neighbor3") {
  ftk::kd<double, 3> mykd;

  std::vector<double> pts(npts * 3);
  for (size_t i = 0; i < npts * 3; i ++)
    pts[i] = (double)rand() / RAND_MAX;

  mykd.set_inputs(npts, &pts[0]);
  mykd.build();

  for (auto batch = 0; batch < npts; batch ++) {
    std::array<double, 3> x = {(double)rand() / RAND_MAX, (double)rand() / RAND_MAX};

    size_t i = mykd.find_nearest(x);
    fprintf(stderr, "x=%f, %f, %f\n", x[0], x[1], x[2]);
    fprintf(stderr, "%zu: %f, %f, %f\n", i, mykd.pts[i][0], mykd.pts[i][1], mykd.pts[i][2]);
    
    size_t j = mykd.find_nearest_naive(x);
    fprintf(stderr, "%zu: %f, %f, %f\n", j, mykd.pts[j][0], mykd.pts[j][1], mykd.pts[j][2]);

    REQUIRE(i == j);
  }
}

#include "main.hh"
