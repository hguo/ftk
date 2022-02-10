#ifndef _FTK_XGC_STREAM_ADIOS2_HH
#define _FTK_XGC_STREAM_ADIOS2_HH

namespace ftk {
using nlohmann::json;

struct xgc_stream_adios2 : public xgc_stream
{
  xgc_stream_adios2(const std::string& path, diy::mpi::communicator comm = MPI_COMM_WORLD) : xgc_stream(path, comm) {}
  
  std::string postfix() const { return ".bp"; }

  std::shared_ptr<ndarray_group> request_step(int step);

  bool read_oneddiag();
  bool advance_timestep();
};

inline bool xgc_stream_adios2::read_oneddiag()
{
  const auto f = oneddiag_filename();
  
  try {
    steps.read_bp(f, "step");
    time.read_bp(f, "time");
    psi_mks.read_bp(f, "psi_mks");

    try {
      e_gc_density_avg.read_bp(f, "e_gc_density_avg");
    } catch (...) {
      warn("cannot read e_gc_density_avg, use _df_1d instead");
      try {
        e_gc_density_avg.read_bp(f, "e_gc_density_df_1d");
      } catch (...) {
        warn("cannot read e_gc_density_avg");
        return false;
      }
    }

    try {
      e_parallel_mean_en_avg.read_bp(f, "e_parallel_mean_en_avg");
    } catch (...) {
      warn("cannot read e_parallel_mean_en_avg, use _df_1d instead");
      try {
        e_parallel_mean_en_avg.read_bp(f, "e_parallel_mean_en_df_1d");
      } catch (...) {
        warn("cannot read e_parallel_mean_en_df_1d");
        return false;
      }
    }

    try {
      e_perp_temperature_avg.read_bp(f, "e_perp_temperature_avg");
    } catch (...) {
      warn("cannot read e_perp_temperature_avg, use _df_1d instead");
      try {
        e_perp_temperature_avg.read_bp(f, "e_perp_temperature_df_1d");
      } catch (...) {
        warn("cannot read e_perp_temperature_df_1d");
        return false;
      }
    }

    Te1d = (e_parallel_mean_en_avg + e_perp_temperature_avg) * (2.0 / 3.0);
    return true;
  } catch (...) {
    warn("cannot read oneddiag file");
    return false;
  }
}

inline std::shared_ptr<ndarray_group> xgc_stream_adios2::request_step(int step)
{
  std::shared_ptr<ndarray_group> g(new ndarray_group);

  const auto f = filename_step(step);

  if (is_directory(f)) {
    try {
      std::shared_ptr<ndarray_base>
        dpot( new ndarray<double> ),
        pot0( new ndarray<double> ),
        potm0( new ndarray<double> ),
        eden( new ndarray<double> );

      dpot->read_bp(f, "dpot");
      pot0->read_bp(f, "pot0");
      potm0->read_bp(f, "potm0");
      eden->read_bp(f, "eden");

      g->set("dpot", dpot);
      g->set("pot0", pot0);
      g->set("potm0", potm0);
      g->set("eden", eden);

      return g;
    } catch (...) {
      warn("error requesting step " + std::to_string(step));
      return NULL;
    }
  } else {
    fatal("cannot open file " + f);
    return NULL;
  }
}

inline bool xgc_stream_adios2::advance_timestep()
{
  std::shared_ptr<ndarray_group> g(new ndarray_group);

  if (current_timestep >= start_timestep + ntimesteps)
    return false;

  const auto current_filename = filename( current_timestep );
  fprintf(stderr, "advancing.., %d, %d, %d, filename=%s\n", start_timestep, current_timestep, ntimesteps, current_filename.c_str());

  if (is_directory( current_filename )) {
    std::shared_ptr<ndarray_base> density( new ndarray<double> );
    density->read_bp(current_filename, "dpot");
    g->set("density", density);
    
    // std::shared_ptr<ndarray_base> Er( new ndarray<double> );
    // Er->read_adios2(current_filename, "Er");
    // g->set("Er", Er);
  
    callback(current_timestep, g);

    current_timestep ++;
    return true;
  } else return false;
}

}

#endif
