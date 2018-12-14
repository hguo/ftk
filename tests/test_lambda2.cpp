#include <ftk/numeric/vortex_criteria.hh>
#include <iostream>

int main(int argc, char **argv)
{
  float J[3][3] = {{0.9102779, 0.44108077, 0.72642273}, 
    {0.39278198, 0.95680469, 0.02683596}, 
    {0.05335823, 0.86960914, 0.43971526}};

  fprintf(stderr, "lambda2=%f\n", ftk::vortex_lambda2_criterion(J));

  return 0;
}
