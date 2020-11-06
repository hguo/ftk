#ifndef _FTK_TRACKER_HH
#define _FTK_TRACKER_HH

#include <ftk/ftk_config.hh>
#include <ftk/filters/filter.hh>
#include <ftk/external/diy/master.hpp>

namespace ftk {

struct tracker : public virtual filter
{
  tracker(diy::mpi::communicator comm) : filter(comm), master(comm) {}
  virtual ~tracker() {}
  
  void set_start_timestep(int t) { start_timestep = t;}
  void set_end_timestep(int t) { end_timestep = t; }
  
  virtual void set_current_timestep(int t) {current_timestep = t;}
  int get_current_timestep() const {return current_timestep;}


protected:
  diy::Master master;
  
protected:
  int start_timestep = 0, 
      end_timestep = std::numeric_limits<int>::max();

  int current_timestep = 0;
};

}

#endif
