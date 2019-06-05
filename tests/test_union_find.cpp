#include <gtest/gtest.h>
#include <ftk/basic/union_find.hh>
#include <ftk/basic/simple_union_find.hh>


class union_find_test : public testing::Test {
public:
  const int nruns = 1;
  const double epsilon = 1e-9;
};

// test (sparse) union-find
TEST_F(union_find_test, union_find) {
  ftk::union_find UF = ftk::union_find(); 
  UF.add("00"); 
  UF.add("11"); 
  UF.add("22"); 
  UF.add("33"); 
  UF.add("111"); 

  UF.unite("11", "22");
  UF.unite("22", "33"); 

  EXPECT_TRUE(UF.same_set("00", "00"));
  EXPECT_TRUE(UF.same_set("11", "22"));
  EXPECT_TRUE(UF.same_set("11", "33"));

  EXPECT_TRUE(!UF.same_set("00", "11"));
  EXPECT_TRUE(!UF.same_set("11", "111"));
}

// test simple union-find
TEST_F(union_find_test, simple_union_find) {
  ftk::simple_union_find UF = ftk::simple_union_find(10); 
  UF.unite(1, 2);
  UF.unite(2, 3); 

  EXPECT_TRUE(UF.same_set(0, 0));
  EXPECT_TRUE(UF.same_set(1, 2));
  EXPECT_TRUE(UF.same_set(1, 3));

  EXPECT_TRUE(!UF.same_set(0, 1));
  EXPECT_TRUE(!UF.same_set(1, 5));
}
