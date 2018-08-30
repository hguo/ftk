#include <ftk/tracking_graph.hh>
#include <fstream>
#include <iostream>

using json = nlohmann::json;

int main(int argc, char **argv)
{
  std::ifstream ifs(argv[1]);
  std::string str;
  str.assign(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());

  json j = json::parse(str);
  // fprintf(stderr, "%s\n", j.dump().c_str());

  ftk::tracking_graph<> tg = j;
  tg.relabel();
  // tr.printComponents();

  // ftk::TransitionLayout::generateDotInputFile(tr, "out.dot");
  // system("dot -Tplain out.dot > out.plain");
  // ftk::TransitionLayout::parseDotOutputFile(tr, "out.plain");

  return 0;
}
