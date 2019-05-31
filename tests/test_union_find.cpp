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
  ftk::sparse_quick_union UF = ftk::sparse_quick_union(); 
  UF.unite("11", "22");
  UF.unite("22", "33"); 

  EXPECT_TRUE(UF.same_set("00", "00"));
  EXPECT_TRUE(UF.same_set("11", "22"));
  EXPECT_TRUE(UF.same_set("11", "33"));

  EXPECT_TRUE(!UF.same_set("00", "11"));
  EXPECT_TRUE(!UF.same_set("11", "111"));
}

// test normal union-find
TEST_F(union_find_test, union_find) {
  ftk::weighted_quick_union UF = ftk::weighted_quick_union(10); 
  UF.unite(1, 2);
  UF.unite(2, 3); 

  EXPECT_TRUE(UF.same_set(0, 0));
  EXPECT_TRUE(UF.same_set(1, 2));
  EXPECT_TRUE(UF.same_set(1, 3));

  EXPECT_TRUE(!UF.same_set(0, 1));
  EXPECT_TRUE(!UF.same_set(1, 5));
}
