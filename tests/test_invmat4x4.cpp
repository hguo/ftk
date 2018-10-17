#include <ftk/numerics/invmat.hh>
#include <iostream>

int main(int argc, char **argv)
{
  float m4x4[] = {0.37116367, 0.16844887, 0.99227088, 0.71275718, 
    0.88786179, 0.19169413, 0.33589513, 0.40073562,
    0.47178621, 0.34809309, 0.46421167, 0.64434074,
    0.9928173, 0.22639279, 0.8466231, 0.48266757};
  float invm4x4[16];

  float det = ftk::invmat4x4(m4x4, invm4x4);

  fprintf(stderr, "det=%f\n", det);
  for (int i=0; i<16; i++) {
    fprintf(stderr, "%f\n", invm4x4[i]);
  }

  return 0;
}

/* expected output: 
det=-0.036813
-0.188140
1.603473
-0.810773
0.028887
-4.662595
-6.674231
5.483859
5.105844
-0.079989
-2.514784
-0.084691
2.319083
2.714262
4.243329
-0.755915
-4.450260
*/
