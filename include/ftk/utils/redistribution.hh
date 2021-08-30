#ifndef _DIYEXT_REDISTRIBUTION_HH
#define _DIYEXT_REDISTRIBUTION_HH

#include <ftk/external/diy/mpi.hpp>
#include <ftk/utils/serialization.hh>
#include <ftk/utils/all2all.hh>
#include <numeric>

namespace diy { namespace mpi {

template <typename K, typename V>
inline void redistribute(const communicator& comm, 
    const std::map<K/*root key*/, std::map<K, V>/*intersections*/>& in,
    std::map<K, std::map<K, V>>& out)
{
  const int np = comm.size(), rank = comm.rank();

  std::vector<std::map<K, std::map<K, V>>> ins(np);
  for (const auto &kv : in) {
    const int rank = (std::hash<K>{}(kv.first)) % np;
    ins[rank].insert(kv);
  }

  // all2all(comm, ins, out);
#if FTK_HAVE_MPI
  // prepare sendbuf
  std::string sendbuf;
  std::vector<int> sendcounts(np), sdispls(np);
  StringBuffer sb(sendbuf);
  for (int i = 0; i < np; i ++) {
    if (false) { // rank == i) {
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
  // fprintf(stderr, "calculating recvcounts...\n");
  std::vector<int> allcounts(np * np);
  MPI_Allgather(&sendcounts[0], np, MPI_INT, 
      &allcounts[0], np, MPI_INT, comm);

  // fprintf(stderr, "counts=%d, %d, %d, %d\n", sendcounts[0], sendcounts[1], sendcounts[2], sendcounts[3]);
  // fprintf(stderr, "sdspls=%d, %d, %d, %d\n", sdispls[0], sdispls[1], sdispls[2], sdispls[3]);
  // fprintf(stderr, "allcounts=%d, %d, %d, %d; %d, %d, %d, %d; %d, %d, %d, %d; %d, %d, %d, %d\n", 
  //    allcounts[0], allcounts[1], allcounts[2], allcounts[3], 
  //    allcounts[4], allcounts[5], allcounts[6], allcounts[7], 
  //    allcounts[8], allcounts[9], allcounts[10], allcounts[11], 
  //    allcounts[12], allcounts[13], allcounts[14], allcounts[15]);

  std::vector<int> recvcounts(np), rdispls(np);
  for (int i = 0; i < np; i ++) {
    rdispls[i] = (i == 0) ? 0 : (rdispls[i-1] + recvcounts[i-1]);
    recvcounts[i] = allcounts[i*np+rank];
  }
  
  // fprintf(stderr, "recvcounts=%d, %d, %d, %d\n", recvcounts[0], recvcounts[1], recvcounts[2], recvcounts[3]);
  // fprintf(stderr, "rdspls=%d, %d, %d, %d\n", rdispls[0], rdispls[1], rdispls[2], rdispls[3]);

  std::string recvbuf;
  recvbuf.resize( std::accumulate(recvcounts.begin(), recvcounts.end(), 0) );

  // all2all
  // fprintf(stderr, "all2all...\n");
  MPI_Alltoallv(&sendbuf[0], &sendcounts[0], &sdispls[0], MPI_CHAR,
      &recvbuf[0], &recvcounts[0], &rdispls[0], MPI_CHAR, comm);
  // fprintf(stderr, "all2all done.\n");

  // unserialize
  // fprintf(stderr, "unserialize... #recvbuf=%zu\n", recvbuf.size());
  // out = ins[rank];
  for (int i = 0; i < np; i ++) {
    // if (i == rank) continue;
    // else {
    StringBuffer bb(recvbuf);
    while (bb) {
      std::map<K, std::map<K, V>> map;
      load(bb, map);
      for (const auto &kv : map)
        out[kv.first].insert(kv.second.begin(), kv.second.end());
    }
  }
#else
  out = in;
#endif
}

#if 0
template <typename K, typename V>
inline void redistribute(const communicator& comm, 
    const std::multimap<K, V>& in, 
    std::multimap<K, V> &out)
{
  const int np = comm.size();

  std::vector<std::multimap<K, V>> ins(np);
  for (const auto &kv : in) {
    const int r = std::hash<K>{}(kv.first) % np;
    ins[r].insert(kv);
  }

  all2all(comm, ins, out);
}
#endif

} // namespace mpi
} // namespace diy

#endif
