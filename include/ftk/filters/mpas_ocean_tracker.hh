#ifndef _FTK_MPAS_OCEAN_TRACKER_HH
#define _FTK_MPAS_OCEAN_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/mesh/simplicial_mpas_2d_mesh.hh>

namespace ftk {

struct mpas_ocean_tracker : public virtual tracker {
  mpas_ocean_tracker(diy::mpi::communicator comm, std::shared_ptr<simplicial_mpas_2d_mesh<>> m_) : 
    tracker(comm), 
    m(m_) {}
  virtual ~mpas_ocean_tracker() {}

public:
  void initialize() {}

protected:
  std::shared_ptr<simplicial_mpas_2d_mesh<>> m;
};

}

#endif
