#ifndef _FTK_XGC_STREAM_H5_HH
#define _FTK_XGC_STREAM_H5_HH

// #include <ftk/io/xgc_stream.hh>

namespace ftk {
using nlohmann::json;

struct xgc_stream_h5 : public xgc_stream
{
  xgc_stream_h5(const std::string& path, diy::mpi::communicator comm = MPI_COMM_WORLD) : xgc_stream(path, comm) {}
  
  std::string postfix() const { return ".h5"; }

  bool advance_timestep();
};

inline bool xgc_stream_h5::advance_timestep()
{
  std::shared_ptr<ndarray_group> g(new ndarray_group);

  if (current_timestep >= start_timestep + ntimesteps)
    return false;

  const auto current_filename = filename( current_timestep );

  if (file_exists( current_filename )) {
    std::shared_ptr<ndarray_base> density( new ndarray<double> );
    density->read_h5(current_filename, "dneOverne0");
    g->set("density", density);
    
    std::shared_ptr<ndarray_base> Er( new ndarray<double> );
    Er->read_h5(current_filename, "Er");
    g->set("Er", Er);
  
    callback(current_timestep, g);

    current_timestep ++;
    return true;
  } else return false;
}


}

#endif
