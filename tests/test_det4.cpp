#include <ftk/numerics/det.hh>
#include <ftk/numerics/print.hh>
#include <iostream>

int main(int argc, char **argv)
{
  float m[] = {
  0.96170121, 0.25145996, 0.15267866, 0.13361515,
  0.72897646, 0.05956134, 0.59597585, 0.00813389,
  0.82357783, 0.50537611, 0.56218735, 0.77733046,
  0.88185366, 0.51499935, 0.1394362,  0.68746057};
  
  ftk::print4x4("m", m);
  fprintf(stderr, "%.12f\n", ftk::det4(m));

  return 0;
}
