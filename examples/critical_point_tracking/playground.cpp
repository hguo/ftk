#include <iostream>
#include <ftk/ndarray.hh>

int main(int argc, char **argv)
{
  ftk::ndarray<double> array;
  array.reshape(128, 128);
  array.copy_to_device();

  return 0;
}
