#define CATCH_CONFIG_RUNNER
#include "catch.hh"
#include <ftk/basic/union_find.hh>
#include <ftk/basic/simple_union_find.hh>
#include <string>

// test (sparse) union-find
TEST_CASE("union_find") {
  ftk::union_find<std::string> UF; 
  UF.add("00"); 
  UF.add("11"); 
  UF.add("22"); 
  UF.add("33"); 
  UF.add("111"); 

  UF.unite("11", "22");
  UF.unite("22", "33"); 

  REQUIRE(UF.same_set("00", "00"));
  REQUIRE(UF.same_set("11", "22"));
  REQUIRE(UF.same_set("11", "33"));

  REQUIRE(!UF.same_set("00", "11"));
  REQUIRE(!UF.same_set("11", "111"));
}

// test simple union-find
TEST_CASE("simple_union_find") {
  ftk::simple_union_find<int> UF(10); 
  UF.unite(1, 2);
  UF.unite(2, 3); 

  REQUIRE(UF.same_set(0, 0));
  REQUIRE(UF.same_set(1, 2));
  REQUIRE(UF.same_set(1, 3));

  REQUIRE(!UF.same_set(0, 1));
  REQUIRE(!UF.same_set(1, 5));
}

int main(int argc, char **argv)
{
  Catch::Session session;
  return session.run(argc, argv);
}
