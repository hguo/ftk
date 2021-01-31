#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include "constants.hh"
#include <ftk/ndarray/stream.hh>
#include <ftk/ndarray/writer.hh>

bool write(const json& jstream, const json& jwriter)
{
  ftk::ndarray_stream<> stream;
  stream.configure(jstream);

  ftk::ndarray_writer<> writer;
  writer.configure(jwriter);
  writer.consume(stream);

  stream.start();
  stream.finish();

  return true; // TODO
}

TEST_CASE("io_write_float64_woven") {
  CHECK(write(js_woven_synthetic, jw_woven_float64));
}

TEST_CASE("io_write_float64_woven_perturbation") {
  CHECK(write(js_woven_synthetic_perturbation, jw_woven_float64_perturbation));
}

TEST_CASE("io_write_float32_tornado") {
  CHECK(write(js_tornado_synthetic, jw_tornado_float32));
}

TEST_CASE("io_write_float32_moving_extremum_3d_synthetic") {
  CHECK(write(js_moving_extremum_3d_synthetic, jw_moving_extremum_3d_float32));
}

#if FTK_HAVE_NETCDF
TEST_CASE("io_write_nc_woven") {
  CHECK(write(js_woven_synthetic, jw_woven_nc_unlimited_time));
}

TEST_CASE("io_write_nc_no_time_woven") {
  CHECK(write(js_woven_synthetic, jw_woven_nc_no_time));
}

//TEST_CASE("io_write_nc_tornado") {
//  CHECK(write(js_tornado_synthetic, jw_tornado_nc));
//}
#endif

#if FTK_HAVE_VTK
TEST_CASE("io_write_vti_woven") {
  CHECK(write(js_woven_synthetic, jw_woven_vti));
}

TEST_CASE("io_write_vti_double_gyre") {
  CHECK(write(js_double_gyre_synthetic, jw_double_gyre_vti));
}

TEST_CASE("io_write_vti_volcano") {
  CHECK(write(js_volcano_2d_synthetic, jw_volcano_vti));
}

TEST_CASE("io_write_vti_moving_extremum_2d_synthetic") {
  CHECK(write(js_moving_extremum_2d_synthetic, jw_moving_extremum_2d_vti));
}

TEST_CASE("io_write_vti_moving_ramp_3d_synthetic") {
  CHECK(write(js_moving_ramp_3d_synthetic, jw_moving_ramp_3d_vti));
}
#endif

#include "main.hh"
