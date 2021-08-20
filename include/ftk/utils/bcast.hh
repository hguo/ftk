#ifndef _DIYEXT_BCAST_HH
#define _DIYEXT_BCAST_HH

#include <ftk/external/diy/mpi.hpp>
#include <ftk/utils/serialization.hh>
#include <numeric>

namespace diy { namespace mpi {

template <typename Obj>
inline void bcastv(const communicator& comm, Obj& inout, int root = 0)
{
#if FTK_HAVE_MPI
  const int np = comm.size(), rank = comm.rank();

  // prepare sendbuf
  int count;
  std::string buf;
  if (comm.rank() == root) {
    serializeToString(inout, buf);
    count = buf.size();
  }

  // bcast size
  MPI_Bcast(&count, 1, MPI_INT, root, comm);
  // fprintf(stderr, "bcast_size=%d\n", count);
  buf.resize(count);

  MPI_Bcast(&buf[0], count, MPI_CHAR, root, comm);

  if (comm.rank() != root) {
    StringBuffer sb(buf);
    load(sb, inout);
  }
#else 
  // nothing to do
#endif
}

}}


#endif
