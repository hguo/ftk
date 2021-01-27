#ifndef _FTK_UNSTRUCTURED_2D_TRACKER_HH
#define _FTK_UNSTRUCTURED_2D_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/mesh/simplicial_unstructured_extruded_2d_mesh.hh>

namespace ftk {

struct unstructured_2d_tracker : public virtual tracker {
  unstructured_2d_tracker(diy::mpi::communicator comm, const simplicial_unstructured_extruded_2d_mesh<> &m) : 
    m(simplicial_unstructured_extruded_2d_mesh<>(m)), tracker(comm) {}
  virtual ~unstructured_2d_tracker() {}

public:
  void initialize() {}

protected:
  const simplicial_unstructured_extruded_2d_mesh<> m;
};

}

#endif
