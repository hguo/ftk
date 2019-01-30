#include <ftk/numeric/interval.hh>
#include <ftk/numeric/interval.hh>
#include <iostream>

int main(int argc, char **argv)
{
  // ftk::basic_interval<int> I(0, 10), I1(9, 30);
  // std::cout << intersect(I, I1) << std::endl;

  ftk::disjoint_intervals<long long int> I(10);
  I.set_to_complete();
  I.join(20, 30);
  I.join(-10, -5);
  I.join(15, 25);
  I.join(30, 60);
  I.join(0, 10);

  std::cout << I << std::endl;

  return 0;
}
