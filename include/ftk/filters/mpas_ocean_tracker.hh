#ifndef _FTK_MPAS_OCEAN_TRACKER_HH
#define _FTK_MPAS_OCEAN_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/mesh/mpas_mesh.hh>

namespace ftk {

struct mpas_ocean_tracker : public virtual tracker {
  mpas_ocean_tracker(diy::mpi::communicator comm, std::shared_ptr<mpas_mesh<>> m_) : 
    tracker(comm), 
    m(m_) {}
  virtual ~mpas_ocean_tracker() {}

public:
  void initialize() {}

protected:
  std::shared_ptr<mpas_mesh<>> m;
};

}

#endif
