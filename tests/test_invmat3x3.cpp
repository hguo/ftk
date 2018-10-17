#include <ftk/numerics/invmat.hh>
#include <iostream>

int main(int argc, char **argv)
{
  float m3x3[] = {0.9102779, 0.44108077, 0.72642273, 0.39278198, 0.95680469, 0.02683596, 0.05335823, 0.86960914, 0.43971526};
  float invm3x3[9];

  float det = ftk::invmat3x3(m3x3, invm3x3);

  fprintf(stderr, "det=%f\n", det);
  for (int i=0; i<9; i++) {
    fprintf(stderr, "%f\n", invm3x3[i]);
  }

  return 0;
}

/* expected output: 
det=0.497218
0.799217
0.880407
-1.374062
-0.344478
0.727051
0.524715
0.584278
-1.544697
1.403228
*/
