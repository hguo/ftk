#include <ftk/numerics/invmat.hh>
#include <iostream>

int main(int argc, char **argv)
{
  float m[] = {0.72758122, 0.31241218, 0.84617905, 0.82793148};
  float inv[4];

  float det = ftk::invmat2(m, inv);

  fprintf(stderr, "det=%f\n", det);
  for (int i=0; i<4; i++) {
    fprintf(stderr, "%f\n", inv[i]);
  }

  return 0;
}

/* expected output: 
det=0.338031
2.449279
-0.924212
-2.503261
2.152411
*/
