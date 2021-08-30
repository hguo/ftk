#ifndef _FTK_CRITICAL_POINT_TRACKER_3D_UNSTRUCTURED_HH
#define _FTK_CRITICAL_POINT_TRACKER_3D_UNSTRUCTURED_HH

#include <ftk/config.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/adjugate.hh>
#include <ftk/numeric/symmetric_matrix.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/geometry/curve2vtk.hh>
#include <ftk/filters/critical_line_tracker.hh>
#include <ftk/filters/unstructured_3d_tracker.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/mesh/simplicial_unstructured_extruded_3d_mesh.hh>
#include <ftk/external/diy/serialization.hpp>

namespace ftk {

struct critical_line_tracker_3d_unstructured :
  public critical_line_tracker, public unstructured_3d_tracker
{
  critical_line_tracker_3d_unstructured(diy::mpi::communicator comm, const simplicial_unstructured_3d_mesh<>& m) : 
    critical_line_tracker(comm), unstructured_3d_tracker(comm, m);
  
  void initialize() {}
  void finalize();
  void reset() {}

  void update_timestep();
;

}
