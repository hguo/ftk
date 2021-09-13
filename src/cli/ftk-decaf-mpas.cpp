#define TRANSPORT_MPI

#include <mpi.h>
#include <ftk/external/diy/mpi.hpp>
#include <decaf/decaf.hpp>
#include <bredala/data_model/simplefield.hpp>

using namespace decaf;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  Workflow workflow;
  Workflow::make_wflow_from_json(workflow, argv[1]);

  Decaf *decaf = new Decaf(MPI_COMM_WORLD, workflow);
  diy::mpi::communicator world(decaf->con_comm_handle());

  fprintf(stderr, "%d\n", world.size());

  delete decaf;

  MPI_Finalize();
  return 0;
}
