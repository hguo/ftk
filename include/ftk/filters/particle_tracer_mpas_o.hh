#ifndef _FTK_PARTICLE_TRACER_MPAS_O_HH
#define _FTK_PARTICLE_TRACER_MPAS_O_HH

#include <ftk/config.hh>
#include <ftk/filters/particle_tracer.hh>

namespace ftk {

struct particle_tracer_mpas_o : public particle_tracer
{
  particle_tracer_mpas_o(diy::mpi::communicator comm) : particle_tracer(comm) {}
  virtual ~particle_tracer_mpas_o() {}

  void update_timestep();
};

//////////

inline void particle_tracer_mpas_o::update_timestep()
{
  if (comm.rank() == 0) fprintf(stderr, "current_timestep=%d\n", current_timestep);
}

}

#endif
