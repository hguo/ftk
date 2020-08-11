#include "constants.hh"
#include <ftk/ndarray.hh>
#include <ftk/mesh/simplicial_unstructured_2d_mesh.hh>

using nlohmann::json;

int main(int argc, char **argv)
{
  diy::mpi::environment env;

#if FTK_HAVE_ADIOS2
  // adios2::ADIOS adios("adios2.xml", MPI_COMM_WORLD);
  adios2::ADIOS adios(MPI_COMM_WORLD);
  adios2::IO io = adios.DeclareIO("BPReader");
  adios2::Engine reader = io.Open(argv[1], adios2::Mode::Read); // , MPI_COMM_SELF);

  const std::map<std::string, adios2::Params> variables = io.AvailableVariables();
  for (const auto pair : variables) {
    std::cerr << pair.first << std::endl;
  }

  ftk::ndarray<double> array;
  array.from_adios2(io, reader, "dpot");

  reader.Close();
#endif

  return 0;
}
