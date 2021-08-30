#ifndef _FTK_MPAS_STREAM_HH
#define _FTK_MPAS_STREAM_HH

#include <ftk/object.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/ndarray_group.hh>

namespace ftk {

using nlohmann::json;

struct mpas_stream : public object {
  mpas_stream(const std::string& path, diy::mpi::communicator comm=MPI_COMM_WORLD);

};

}

#endif
