#include <ftk/numeric/interval.hh>
#include <ftk/numeric/interval.hh>
#include <iostream>

int main(int argc, char **argv)
{
  ftk::basic_interval<float> my_interval(0, 10);

  std::cout << my_interval.contains(3) << std::endl;
  std::cout << my_interval.contains(10) << std::endl;
  std::cout << my_interval.contains(18) << std::endl;
 
#if 0
  ftk::basic_interval<float> my_interval1(20, true, 30, false);
  std::cout << my_interval1.contains(30) << std::endl;
  std::cout << my_interval1.empty() << std::endl;
#endif
  
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

  return 0;
}
