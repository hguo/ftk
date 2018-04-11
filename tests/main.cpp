#include <ftk/transition/transition.h>
#include <ftk/transition/layout.h>
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

  ftk::Transition tr = j;
  // tr.relabel();
  // tr.printComponents();

  ftk::TransitionLayout::generateDotFile(tr, argv[2]);

  return 0;
}
