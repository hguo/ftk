#ifndef _FTK_XGC_POINCARE_FILTER_HH
#define _FTK_XGC_POINCARE_FILTER_HH

#include <ftk/config.hh>
#include <ftk/filters/xgc_tracker.hh>
#include <ftk/features/feature_line.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/clamp.hh>

#if FTK_HAVE_CUDA // CUDA is actually required for this filter
#include "xgc_blob_filament_tracker.cuh"
#endif

namespace ftk {

struct xgc_poincare_filter : public filter
{
  xgc_poincare_filter(diy::mpi::communicator comm, 
      std::shared_ptr<simplicial_xgc_3d_mesh<>> mx);
  
  void set_use_deltaB(bool);
  void set_revolutions(int);

  void initialize_regular_seeds(
      int nseeds_per_sector = 5000,
      int nsectors = 1,
      double psin0 = 0.8, 
      double psin1 = 1.03);

protected:
  bool use_deltaB = true;
  int revolutions = 3000;
};

}

#endif
