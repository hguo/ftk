#include <ftk/numerics/det.hh>
#include <iostream>

int main(int argc, char **argv)
{
  float m[] = {0.9102779, 0.44108077, 0.72642273, 0.39278198, 0.95680469, 0.02683596, 0.05335823, 0.86960914, 0.43971526};
  fprintf(stderr, "%f\n", ftk::det3(m));

  return 0;
}

/* expected output:
0.497218
*/
