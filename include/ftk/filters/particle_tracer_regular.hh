#ifndef _FTK_PARTICLE_TRACER_REGULAR_HH
#define _FTK_PARTICLE_TRACER_REGULAR_HH

#include <ftk/ndarray.hh>
#include <ftk/mesh/lattice_partitioner.hh>
#include <ftk/filters/particle_tracer.hh>
#include <ftk/filters/regular_tracker.hh>
#include <ftk/utils/gather.hh>

namespace ftk {

struct particle_tracer_regular : public particle_tracer, public regular_tracker
{
  particle_tracer_regular(diy::mpi::communicator comm, int nd) : 
    particle_tracer(comm), 
    regular_tracker(comm, nd), 
    tracker(comm) 
  {} 

  virtual ~particle_tracer_regular() {}
};

}

#endif
