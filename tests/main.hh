#include <ftk/external/diy/mpi.hpp>

int main(int argc, char **argv)
{
  int requested = MPI_THREAD_FUNNELED, provided;
#if FTK_HAVE_MPI
  MPI_Init_thread(&argc, &argv, requested, &provided);
#endif

  Catch::Session session;
  int result = session.run(argc, argv);

#if FTK_HAVE_MPI
  MPI_Finalize();
#endif

  return result;
}
