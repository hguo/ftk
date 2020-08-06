#include "constants.hh"
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>

using nlohmann::json;

int main(int argc, char **argv)
{
  diy::mpi::environment env;
  
  ftk::simplicial_unstructured_2d_mesh<> mesh;
  mesh.from_vtk_file(argv[1]);
}
