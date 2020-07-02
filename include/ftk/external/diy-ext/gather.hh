#ifndef _DIYEXT_GATHER_HH
#define _DIYEXT_GATHER_HH

#include <ftk/external/diy/mpi.hpp>
#include <ftk/external/diy-ext/serialization.hh>

namespace diy { namespace mpi {

template <typename K, typename V> // merging an std::map to root using gather
inline void gather(const communicator& comm, const std::map<K, V>& in, std::map<K, V> &out, int root)
{
#if FTK_HAVE_MPI
  // serialize input map
  std::string serialized_in;
  if (comm.rank() != root) // avoid serializing data from the root proc
    serializeToString(in, serialized_in);

  // gathering length of serialized data
  int length_serialized_in = serialized_in.size();
  std::vector<int> all_length_serialized_in(comm.size(), 0);
  MPI_Gather(&length_serialized_in, 1, MPI_INT, 
      &all_length_serialized_in[0], 1, MPI_INT, 
      root, comm);

  // preparing buffer and displacements for MPI_Gatherv
  std::string buffer;
  std::vector<int> displs;
  if (comm.rank() == root) {
    buffer.resize(std::accumulate( // prepare buffer
          all_length_serialized_in.begin(), 
          all_length_serialized_in.end(), 0));

    displs.resize(all_length_serialized_in.size(), 0); // prepare displs
    for (int i = 1; i < comm.size(); i ++)
      displs[i] = displs[i-1] + all_length_serialized_in[i-1];
  }

  // call MPI_Gatherv
  MPI_Gatherv(serialized_in.data(), serialized_in.size(), MPI_CHAR, 
      &buffer[0], all_length_serialized_in.data(), displs.data(), MPI_CHAR, 
      root, comm);

  // unserailization (root proc only)
  if (comm.rank() == root) {
    out = in;
    StringBuffer sb(buffer);
    while (sb) {
      std::map<K, V> map;
      load(sb, map);
      for (const auto &kv : map)
        out.insert(kv);
    }
  }
#else
  out = in;
#endif
}

}
}

#endif
