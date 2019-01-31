#include <ftk/numeric/disjoint_intervals.hh>
#include <iostream>

int main(int argc, char **argv)
{
#if 0
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
#endif


#if 1
  ftk::disjoint_intervals<float> I(-200, -100);
  I.join(100, 200);
  I.join(ftk::basic_interval<float>::lower_inf(), -300);
  std::cerr << "I=" << I << std::endl;

  ftk::disjoint_intervals<float> J(-400, -150);
  J.join(150, ftk::basic_interval<float>::upper_inf());
  std::cerr << "J=" << J << std::endl;

  J.intersect(I);

  std::cerr << "intersect(I,J)=" << J << std::endl;
#endif

#if 0
  ftk::basic_interval<float> i(ftk::basic_interval<float>::lower_inf(), -300);
  ftk::basic_interval<float> j(150, ftk::basic_interval<float>::upper_inf());
  // std::cerr << "intersect(i,j)=" << intersect(i, j) << std::endl;
  i.intersect(j);
  std::cerr << i << std::endl;
#endif


#if 0
  ftk::disjoint_intervals<float> K;
  K.set_to_complete();
  K.intersect(I);
  std::cerr << I2 << std::endl;
#endif

  return 0;
}
