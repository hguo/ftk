#ifndef _DIYEXT_SCATTER_HH
#define _DIYEXT_SCATTER_HH

#include <ftk/external/diy/mpi.hpp>
#include <ftk/utils/serialization.hh>
#include <numeric>

namespace diy { namespace mpi {

template <typename Obj>
inline void scatterv(const communicator& comm, 
    const std::vector<Obj> &ins,
    Obj &out, 
    int root = 0)
{
#if FTK_HAVE_MPI
  assert(ins.size() == comm.size());

  const int np = comm.size(), rank = comm.rank();

  // prepare sendbuf
  std::string sendbuf;
  StringBuffer sb(sendbuf);
  std::vector<int> sendcounts(np), displs(np);
  for (int i = 0; i < ins.size(); i ++) {
    displs[i] = sendbuf.size();
    save(sb, ins[i]);
    sendcounts[i] = sendbuf.size() - displs[i];
  }

  // bcast counts and displs
  broadcast(comm, sendcounts, root);
  broadcast(comm, displs, root);

  // prepare recvbuf
  std::string recvbuf;
  recvbuf.resize( sendcounts[rank] );

  // call MPI_Scatterv
  MPI_Scatterv(&sendbuf[0], &sendcounts[0], &displs[0], MPI_CHAR,
      &recvbuf[0], recvbuf.size(), MPI_CHAR, 
      root, comm);

  // unseralize
  StringBuffer bb(recvbuf);
  load(bb, out);
#else
  out = ins[0];
#endif
}

}}


#endif
