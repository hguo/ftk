#ifndef _FTK_XGC_STREAM_ADIOS2_HH
#define _FTK_XGC_STREAM_ADIOS2_HH

namespace ftk {
using nlohmann::json;

struct xgc_stream_adios2 : public xgc_stream
{
  xgc_stream_adios2(const std::string& path, diy::mpi::communicator comm = MPI_COMM_WORLD) : xgc_stream(path, comm) {}
  
  std::string postfix() const { return ".bp"; }

  std::shared_ptr<ndarray_group> request_step(int step);

  bool read_units();

  bool read_oneddiag();
  bool advance_timestep();
};

inline bool xgc_stream_adios2::read_units()
{
  xgc_units_t u;
  bool succ = true;

  if (!u.read(units_filename())) 
    succ = u.read_bp( path + "/xgc.units.bp");

  if (succ) 
    m2->set_units(u);

  // iphi = std::max(1, u.sml_wedge_n);
  return succ;
}

inline bool xgc_stream_adios2::read_oneddiag()
{
  const auto f = oneddiag_filename();
  
  try {
    steps.read_bp(f, "step", NDARRAY_ADIOS2_STEPS_ALL);
    for (int i = 0; i < steps.size(); i ++) {
      step2istep[steps[i]] = i;
      // fprintf(stderr, "step=%d, istep=%d\n", steps[i], i);
    }

    time.read_bp(f, "time");

#if 0
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
#endif
    return true;
  } catch (...) {
    warn("cannot read oneddiag file");
    return false;
  }
}

inline std::shared_ptr<ndarray_group> xgc_stream_adios2::request_step(int step)
{
  std::shared_ptr<ndarray_group> g(new ndarray_group);

  const int istep = step2istep[step];
  fprintf(stderr, "istep=%d\n", istep);

  // read f1d
  const auto f1d = oneddiag_filename();
  const auto f = filename_step(step);
  if (is_directory(f1d) && is_directory(f)) {
    try {
      std::shared_ptr<ndarray_base> // 1d profiles
        psi_mks( new ndarray<double> ),
        e_gc_density_avg( new ndarray<double> ), 
        e_parallel_mean_en_avg( new ndarray<double> ),
        e_perp_temperature_avg( new ndarray<double> );

      psi_mks->read_bp(f1d, "psi_mks", istep);
      g->set("psi_mks", psi_mks);
    
      try {
        e_gc_density_avg->read_bp(f1d, "e_gc_density_avg", istep);
      } catch (...) {
        warn("cannot read e_gc_density_avg, use _df_1d instead");
        try {
          e_gc_density_avg->read_bp(f1d, "e_gc_density_df_1d", istep);
        } catch (...) {
          warn("cannot read e_gc_density_df_1d");
        }
      }
      g->set("e_gc_density_avg", e_gc_density_avg);

      try {
        e_parallel_mean_en_avg->read_bp(f1d, "e_parallel_mean_en_avg", istep);
      } catch (...) {
        warn("cannot read e_parallel_mean_en_avg, use _df_1d instead");
        try {
          e_parallel_mean_en_avg->read_bp(f1d, "e_parallel_mean_en_df_1d", istep);
        } catch (...) {
          warn("cannot read e_parallel_mean_en_df_1d");
        }
      }
      g->set("e_parallel_mean_en_avg", e_parallel_mean_en_avg);

      try {
        e_perp_temperature_avg->read_bp(f1d, "e_perp_temperature_avg", istep);
      } catch (...) {
        warn("cannot read e_perp_temperature_avg, use _df_1d instead");
        try {
          e_perp_temperature_avg->read_bp(f1d, "e_perp_temperature_df_1d", istep);
        } catch (...) {
          warn("cannot read e_perp_temperature_df_1d");
        }
      }
      g->set("e_perp_temperature_avg", e_perp_temperature_avg);

      std::shared_ptr<ndarray_base> // 3d
        dpot( new ndarray<double> ),
        pot0( new ndarray<double> ),
        potm0( new ndarray<double> ),
        eden( new ndarray<double> );

      dpot->read_bp(f, "dpot");
      g->set("dpot", dpot);
      
      pot0->read_bp(f, "pot0");
      g->set("pot0", pot0);
      
      potm0->read_bp(f, "potm0");
      g->set("potm0", potm0);
      
      eden->read_bp(f, "eden");
      g->set("eden", eden);

      ndarray<double> dneOverne0 = mx3->derive_turbulence(
          g->get<double>("dpot"),
          g->get<double>("pot0"), 
          g->get<double>("potm0"),
          g->get<double>("eden"),
          g->get<double>("psi_mks"),
          g->get<double>("e_gc_density_avg"),
          g->get<double>("e_perp_temperature_avg"),
          g->get<double>("e_parallel_mean_en_avg"));
      g->set("dneOverne0", dneOverne0);

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
