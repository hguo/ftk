#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <vector>
#include <cmath>
#include <ftk/ndarray.hh>

const int nruns = 100000;
const double epsilon = 1e-9;

TEST_CASE("2D_lerp_2") 
{
  ftk::ndarray<double> A({2, 2, 2}); // with 2 components
  A.from_vector({0.0, 1.0, 
                 1.0, 2.0, 
                 2.0, 3.0, 
                 3.0, 4.0});
  A.set_multicomponents(1);

  double x[2] = {0.5, 0.5};
  double v[2];

  A.mlerp(x, v);
    
  REQUIRE(v[0] == Approx(1.5).margin(epsilon));
  REQUIRE(v[1] == Approx(2.5).margin(epsilon));
}

TEST_CASE("2D_lerp_1") 
{
  ftk::ndarray<double> A({2, 2}); 
  A.from_vector({0.0, 1.0, 2.0, 3.0});

  double x[2] = {0.5, 0.5};
  double v[1];

  A.mlerp(x, v);
    
  REQUIRE(v[0] == Approx(1.5).margin(epsilon));
}

#include "main.hh"
