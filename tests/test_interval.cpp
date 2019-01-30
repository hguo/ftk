#include <ftk/numeric/interval.hh>
#include <ftk/numeric/interval.hh>
#include <iostream>

int main(int argc, char **argv)
{
  ftk::basic_interval<int> I(0, 10), I1(9, 30);

  std::cout << intersect(I, I1) << std::endl;
#if 0
  std::cout << my_interval.contains(3) << std::endl;
  std::cout << my_interval.contains(10) << std::endl;
  std::cout << my_interval.contains(18) << std::endl;
  
  ftk::basic_interval<float> my_interval2(1); 
  std::cout << my_interval2.empty() << std::endl;
  std::cout << my_interval2.singleton() << std::endl;

  ftk::basic_interval<float> my_interval3;
  my_interval3.set_to_empty();
  std::cout << my_interval3.empty() << std::endl;
  my_interval3.set_to_complete();
  std::cout << my_interval3.empty() << std::endl;
  std::cout << my_interval3.contains(1000) << std::endl;
  std::cout << my_interval3.complete() << std::endl;
#endif
  return 0;
}
