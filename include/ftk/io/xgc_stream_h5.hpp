#ifndef _FTK_XGC_STREAM_H5_HH
#define _FTK_XGC_STREAM_H5_HH

// #include <ftk/io/xgc_stream.hh>

namespace ftk {
using nlohmann::json;

struct xgc_stream_h5 : public xgc_stream
{
  xgc_stream_h5(const std::string& path, diy::mpi::communicator comm = MPI_COMM_WORLD) : xgc_stream(path, comm) {}
  
  std::string postfix() const { return ".h5"; }
  
  std::shared_ptr<ndarray_group> request_step(int step) { return NULL; } // TODO

  bool read_oneddiag();
  bool advance_timestep();
};

inline bool xgc_stream_h5::read_oneddiag()
{
  const auto f = oneddiag_filename();
  ndarray<double> etemp_par, etemp_per;
  
  steps.read_h5(f, "steps");
  time.read_h5(f, "time");
  etemp_par.read_h5(f, "e_parallel_mean_en_avg");
  etemp_per.read_h5(f, "e_perp_temperature_avg");

  Te1d = (etemp_par + etemp_per) * (2.0 / 3.0);

  if (steps.size() == 0) return false;
  else return true;
}

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
