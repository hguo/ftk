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

#if FTK_HAVE_NETCDF
TEST_CASE("io_write_nc_woven") {
  CHECK(write(js_woven_synthetic, jw_woven_nc_unlimited_time));
}

TEST_CASE("io_write_nc_no_time_woven") {
  CHECK(write(js_woven_synthetic, jw_woven_nc_no_time));
}
#endif

#if FTK_HAVE_VTK
TEST_CASE("io_write_vti_woven") {
  CHECK(write(js_woven_synthetic, jw_woven_vti));
}
#endif

int main(int argc, char **argv)
{
  Catch::Session session;
  return session.run(argc, argv);
}
