#ifndef _FTK_UNSTRUCTURED_3D_TRACKER_HH
#define _FTK_UNSTRUCTURED_3D_TRACKER_HH

#include <ftk/config.hh>
#include <ftk/filters/tracker.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/mesh/simplicial_unstructured_extruded_3d_mesh.hh>

namespace ftk {

struct unstructured_3d_tracker : public virtual tracker {
  unstructured_3d_tracker(diy::mpi::communicator comm, std::shared_ptr<simplicial_unstructured_3d_mesh<>> m3_) : 
    m3(m3_),
    m(new simplicial_unstructured_extruded_3d_mesh<>(*m3_)), tracker(comm) {}
  virtual ~unstructured_3d_tracker() {}

public:
  void initialize() {}

protected:
  std::shared_ptr<simplicial_unstructured_3d_mesh<>> m3;
  std::shared_ptr<simplicial_unstructured_extruded_3d_mesh<>> m;
};

}

#endif
