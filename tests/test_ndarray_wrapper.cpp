#include <ftk/ndarray_wrapper.hh>

using namespace ftk;

int main(int argc, char **argv)
{
  ftk::ndarray_wrapper w(NDARRAY_TYPE_DOUBLE);
  std::cerr << w.nd() << std::endl;
}
