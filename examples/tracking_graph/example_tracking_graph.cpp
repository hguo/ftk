#include <ftk/tracking_graph.hh>
#include <fstream>
#include <iostream>

using json = nlohmann::json;

int main(int argc, char **argv)
{
  ftk::tracking_graph<> tg;

  // timestep 0
  tg.add_node(0, 0);
  tg.add_node(0, 1);
  tg.add_node(0, 2);
  
  // timestep 1
  tg.add_node(1, 0);
  tg.add_node(1, 1);

  // timestep 2
  tg.add_node(2, 0);
  tg.add_node(2, 1);
  tg.add_node(2, 2);

  // interval (0, 1)
  tg.add_edge(0, 0, 1, 0);
  tg.add_edge(0, 1, 1, 0);
  tg.add_edge(0, 2, 1, 1);

  // interval (0, 2)
  tg.add_edge(1, 0, 2, 0);
  tg.add_edge(1, 1, 2, 1);
  tg.add_edge(1, 1, 2, 2);

  // update global labels
  tg.relabel();

  // serialize the tracking graph to json (WIP)
  // json j = tg;
  // std::string str = j.dump();
  // std::cout << str << std::endl;

  // parse the tracking graph from json (WIP)
  // json j1 = json::parse(str);
  // ftk::tracking_graph<> tg1 = j1;

  // store the tracking graph in the dot format
  tg.generate_dot_file("dot");
  std::cout << "run the following command to generate the tracking graph visualization in SVG format:" << std::endl
            << "$ dot -Tsvg dot > dot.svg" << std::endl;

  return 0;
}
