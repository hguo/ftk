#include <iostream>
#include <ftk/ndarray/stream.hh>

using nlohmann::json;

int main(int argc, char **argv)
{
  ftk::ndarray_stream<> stream;
  stream.set_input_source_json_file(argv[1]);
  std::cerr << stream.get_json() << std::endl;

  stream.set_callback([&](int k, ftk::ndarray<double> array) {
    fprintf(stderr, "got timestep k=%d\n", k);
    // std::cerr << array << std::endl;
  });

  stream.start();
  stream.finish();

  return 0;
}
