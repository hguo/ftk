#include <gtest/gtest.h>
#include <ftk/basic/union_find.hh>
#include <ftk/basic/sparse_union_find.hh>


class union_find_test : public testing::Test {
public:
  const int nruns = 1;
  const double epsilon = 1e-9;
};

// test sparse union-find
TEST_F(union_find_test, sparse_union_find) {
  ftk::sparse_weighted_quick_union UF = ftk::sparse_weighted_quick_union(); 
  UF.unite(1, 2);
  UF.unite(2, 3); 

  size_t id_0 = UF.root(0); 
  size_t id_1 = UF.root(1); 
  size_t id_2 = UF.root(2); 
  size_t id_3 = UF.root(3); 

  EXPECT_NEAR(0.0, id_0 - 0, epsilon);
  EXPECT_NEAR(0.0, id_1 - id_2, epsilon);
  EXPECT_NEAR(0.0, id_1 - id_3, epsilon);
}

// test normal union-find
TEST_F(union_find_test, union_find) {
  ftk::weighted_quick_union UF = ftk::weighted_quick_union(10); 
  UF.unite(1, 2);
  UF.unite(2, 3); 

  size_t id_0 = UF.root(0); 
  size_t id_1 = UF.root(1); 
  size_t id_2 = UF.root(2); 
  size_t id_3 = UF.root(3); 

  EXPECT_NEAR(0.0, id_0 - 0, epsilon);
  EXPECT_NEAR(0.0, id_1 - id_2, epsilon);
  EXPECT_NEAR(0.0, id_1 - id_3, epsilon);
}
