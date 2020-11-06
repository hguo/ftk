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
  
  virtual int cpdims() const = 0; // featutre dimension
  
  void set_start_timestep(int t) { start_timestep = t;}
  void set_end_timestep(int t) { end_timestep = t; }
  
  virtual void set_current_timestep(int t) {current_timestep = t;}
  int get_current_timestep() const {return current_timestep;}
 
  void set_input_array_partial(bool b) {is_input_array_partial = b;}
  void set_use_default_domain_partition(bool b) {use_default_domain_partition = true;}
 
public:
  virtual void initialize() = 0;
  virtual void finalize() = 0;
  
  virtual bool advance_timestep() = 0;
  virtual void update_timestep() = 0;

protected:
  diy::Master master;
  
protected:
  int start_timestep = 0, 
      end_timestep = std::numeric_limits<int>::max();

  int current_timestep = 0;
 
protected:
  bool is_input_array_partial = false;
  bool use_default_domain_partition = true;
};

}

#endif
