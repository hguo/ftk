#include "constants.hh"
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>

using nlohmann::json;

int main(int argc, char **argv)
{
  diy::mpi::environment env;
  
  ftk::simplicial_unstructured_2d_mesh<> mesh;
  // mesh.from_legacy_vtk_file(argv[1]);
  // mesh.to_vtk_unstructured_grid_file(argv[2]);
  mesh.from_vtk_unstructured_grid_file(argv[1]);

  return 0;
}
