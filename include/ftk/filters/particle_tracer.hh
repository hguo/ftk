#ifndef _FTK_PARTICLE_TRACER_HH
#define _FTK_PARTICLE_TRACER_HH

#include <ftk/config.hh>
#include <ftk/filters/tracker.hh>

namespace ftk {

struct particle_tracer : virtual public tracker
{
  particle_tracer(diy::mpi::communicator comm) : tracker(comm) {}
  virtual ~particle_tracer() {}

  virtual void set_particles(const ndarray<double> &p) { particles = p; }

protected:
  ndarray<double> particles;
};

//// 
}

#endif
