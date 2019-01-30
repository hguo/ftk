#include <ftk/numeric/disjoint_intervals.hh>
#include <iostream>

int main(int argc, char **argv)
{
  // ftk::basic_interval<int> I(0, 10), I1(9, 30);
  // std::cout << intersect(I, I1) << std::endl;

  ftk::disjoint_intervals<float> I(10);
  // ftk::disjoint_intervals<float> I(10);
  // I.set_to_complete();
#if 1
  I.join(20, 30);
  I.join(-10, -5);
  I.join(15, 25);
  I.join(30, 60);
  I.join(0, 10);
  I.join(ftk::basic_interval<float>::lower_inf(), -100);
#endif
  I.join(100, ftk::basic_interval<float>::upper_inf());
  // I.intersect(15);

  // I.intersect(ftk::basic_interval<long long int>::complete_interval());

  std::cout << I << std::endl;

  return 0;
}
