#ifndef _FTK_UNSTRUCTURED_2D_TRACKER_HH
#define _FTK_UNSTRUCTURED_2D_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/mesh/simplicial_unstructured_extruded_2d_mesh.hh>

namespace ftk {

struct unstructured_2d_tracker : public virtual tracker {
  unstructured_2d_tracker(diy::mpi::communicator comm, 
      std::shared_ptr<simplicial_unstructured_extruded_2d_mesh<>> m_) : 
    m(m_), tracker(comm) {}

  virtual ~unstructured_2d_tracker() {}

public:
  void initialize() {}

protected:
  std::shared_ptr<simplicial_unstructured_extruded_2d_mesh<>> m;
};

}

#endif
