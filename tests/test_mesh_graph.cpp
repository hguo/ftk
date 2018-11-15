#include <ftk/mesh_graph/mesh_graph_regular_3d.hh>
#include <ftk/mesh_graph/mesh_graph_regular_3d_tets.hh>

int main(int argc, char **argv)
{
  ftk::mesh_graph_regular_3d_tets<> mg(128, 128, 128);

  fprintf(stderr, "%lu\n", mg.n(2));

  return 0;
}
