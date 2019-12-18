#include <iostream>
#include <ftk/ndarray.hh>

int main(int argc, char **argv)
{
  const int m = 7, n = 8; 
  ftk::ndarray<int> array({m, n}); // n rows and m columns
  for (int i = 0; i < array.nelem(); i ++)
    array[i] = i;
  
  for (int j = 0; j < n; j ++) {
    for (int i = 0; i < m; i ++) 
      fprintf(stderr, "%d, ", array(i, j));
    fprintf(stderr, "\n");
  }

  const int p = 3, q = 2;
  auto subarray = array.slice({3, 2}, {p, q});
  
  for (int j = 0; j < q; j ++) {
    for (int i = 0; i < p; i ++) 
      fprintf(stderr, "%d, ", subarray(i, j));
    fprintf(stderr, "\n");
  }

  return 0;
}
