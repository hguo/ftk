#include <ftk/basic/union_find.h>
#include <ftk/basic/contour_tree.h>
#include <ftk/mesh_graph/regular_2d.h>
#include <ftk/algorithms/sweep_and_merge.h>

int main(int argc, char **argv)
{
#if 1 
  const int W = 5, H = 3;
  std::vector<double> values = {
    0, 15, 4, 7, 2, 
    17, 20, 5, 10, 8, 
    13, 16, 1, 6, 3};
#endif
#if 0
  const int W = 3, H = 3;
  std::vector<double> values = {
    8, 7, 4, 
    6, 5, 2, 
    3, 1, 0};
#endif

  auto ct = ftk::build_contour_tree<size_t, double>(W*H, values, 
      // std::bind(ftk::Get6Neighbors2DRegular<size_t>, W, H, std::placeholders::_1));
      std::bind(&ftk::regular_mesh_2d<size_t>::get6neighbors, W, H, std::placeholders::_1));

  // ct.print_with_values<double>( [&values](size_t i) {return values[i];} );
  ct.reduce();
  ct.print();

  // fprintf(stderr, "====\n");
  // ct.reduce();
  // ct.print_with_values<double>( [&values](size_t i) {return values[i];} );

  return 0;
}
