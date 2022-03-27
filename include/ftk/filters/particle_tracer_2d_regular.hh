#ifndef _FTK_PARTICLE_TRACER_2D_REGULAR_HH
#define _FTK_PARTICLE_TRACER_2D_REGULAR_HH

#include <ftk/config.hh>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/clamp.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
// #include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/adjugate.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/geometry/curve2vtk.hh>
#include <ftk/filters/particle_tracer_regular.hh>
#include <ftk/ndarray.hh>
#include <ftk/ndarray/grad.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/utils/gather.hh>

namespace ftk {

struct particle_tracer_2d_regular : public particle_tracer_regular
{
  particle_tracer_2d_regular(diy::mpi::communicator comm) : 
    particle_tracer_regular(comm, 2), 
    tracker(comm)
  {}
  virtual ~particle_tracer_2d_regular() {}
};

}

#endif
