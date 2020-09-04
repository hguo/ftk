#include "constants.hh"
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>
#include <iostream>

using nlohmann::json;

int main(int argc, char **argv)
{
  diy::mpi::environment env;
  
  ftk::simplicial_unstructured_2d_mesh<> mesh;
  // mesh.from_legacy_vtk_file(argv[1]);
  // mesh.to_vtk_unstructured_grid_file(argv[2]);
  mesh.from_vtk_unstructured_grid_file(argv[1]);

  for (int i = 0; i < 32; i ++) {
    char filename[1024];
    snprintf(filename, 1024, "woven-%03d.vtu\n", i);

    auto woven = ftk::synthetic_woven_2D_unstructured<double>(mesh.get_coords(), i*0.1);
    mesh.scalar_to_vtk_unstructured_grid_data_file(filename, "scalar", woven);
  }

  return 0;
}
