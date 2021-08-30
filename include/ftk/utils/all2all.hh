#ifndef _DIYEXT_ALL2ALL_HH
#define _DIYEXT_ALL2ALL_HH

#include <ftk/external/diy/mpi.hpp>
#include <ftk/utils/serialization.hh>
#include <numeric>

namespace diy { namespace mpi {

template <typename Map> // map = std::map, std::set, or other containers
inline void all2all(const communicator& comm, const std::vector<Map>& ins, Map &out)
{
#if FTK_HAVE_MPI
  const int np = comm.size(), rank = comm.rank();

  // prepare sendbuf
  std::string sendbuf;
  std::vector<int> sendcounts(np), sdispls(np);
  StringBuffer sb(sendbuf);
  for (int i = 0; i < np; i ++) {
    if (rank == i) {
      // out = ins[i];
      // out.insert(ins[i].begin(), ins[i].end());
      sdispls[i] = sendbuf.size();
      sendcounts[i] = 0;
    } else {
      sdispls[i] = sendbuf.size();
      save(sb, ins[i]);
      sendcounts[i] = sendbuf.size() - sdispls[i];
    }
  }

  // calculate recvcounts
  std::vector<int> allcounts(np * np);
  MPI_Allgather(&sendcounts[0], np, MPI_INT, 
      &allcounts[0], np, MPI_INT, comm);

  std::vector<int> recvcounts(np), rdispls(np);
  for (int i = 0; i < np; i ++) {
    rdispls[i] = recvcounts[i];
    recvcounts[i] = allcounts[i*np+rank];
  }

  std::string recvbuf;
  recvbuf.resize( std::accumulate(recvcounts.begin(), recvcounts.end(), 0) );

  // all2all
  MPI_Alltoallv(&sendbuf[0], &sendcounts[0], &sdispls[0], MPI_CHAR,
      &recvbuf[0], &recvcounts[0], &rdispls[0], MPI_CHAR, comm);

  // unserialize
  out = ins[rank];
  StringBuffer bb(recvbuf);
  while (bb) {
    Map map;
    load(bb, map);
    for (const auto &kv : map)
      out.insert(kv);
  }
#endif
}

} // namespace mpi
} // namespace diy

#endif
