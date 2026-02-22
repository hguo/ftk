#ifndef _DIYEXT_BCAST_HH
#define _DIYEXT_BCAST_HH

#include <ftk/external/diy/mpi.hpp>
#include <yaml-cpp/yaml.h>
#include <ftk/utils/serialization.hh>
#include <numeric>

namespace diy { namespace mpi {

using json = YAML::Node;

inline void bcastj(const communicator& comm, json& j, int root = 0)
{
  if (comm.size() == 1) return;

  std::string yaml_str;

  if (comm.rank() == root) {
    YAML::Emitter emitter;
    emitter << j;
    yaml_str = emitter.c_str();
  }

  diy::mpi::broadcast(comm, yaml_str, root);

  if (comm.rank() != root)
    j = YAML::Load(yaml_str);
}

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
