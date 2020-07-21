#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <ftk/ndarray/stream.hh>
#include <ftk/ndarray/writer.hh>

using nlohmann::json;

const json js_woven = {
  {"type", "synthetic"},
  {"name", "woven"}
};

const json js_woven_bin = {
  {"type", "file"},
  {"format", "float64"},
  {"filename", "woven-%04d.bin"}
};

const json jw_woven = {
    {"nd", 2},
    {"format", "float64"},
    {"filename", "woven-%04d.bin"},
    {"variable", "scalar"}
  };

TEST_CASE("io_write_float64_woven_read_float64_woven") {
  ftk::ndarray_stream<> stream;
  stream.configure(js_woven);

  ftk::ndarray_writer<> writer;
  writer.configure(jw_woven);
  writer.consume(stream);

  stream.start();
  stream.finish();
}

int main(int argc, char **argv)
{
  Catch::Session session;
  return session.run(argc, argv);
}
